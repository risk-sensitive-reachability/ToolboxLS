function [ data, g, data0 ] = animateSpinStar(filename, accuracy, compress)
% animateSpinStar: create an animation of the growth of the spinning star.
%
%   [ data, g, data0 ] = animateSpinStar(filename, accuracy, compress)
%  
% This file generates an animation showing how the spinning star grows
%   as time progresses.  It is basically a combination of the files:
%        spinStarDemo (which sets up and handles the level set calculation),
%        spinAnimation (which has the code necessary for animations).
%
% NOTES:
%
% 1) This function will probably only work in Windows, since it uses
%    the avifile command and the avi animation format.
%
% 2) While the animation is being generated, don't move anything
%    (mouse, other windows) in front of the matlab figure window,
%    otherwise that image will be captured into the avi.
%
% 3) If this file stops because of an error, the avi file will remain open.
%    In that case, you must issue a "clear all" command to get rid of
%    the open file handle and then delete the partially completed avi file.
%
%  
% Recreates figure 6.2 from O&F chapter 6, showing a combination of motion 
%   by surface normal and convective rotation of a star-shaped interface 
%   in two dimensions.
%
% The caption claim of "rigid body rotation" is clearly incorrect, since
%   the tips of the star are rotating faster than the center.
%
% This function allows either recreation of the figure or true rigid
%   body rotation (which is not quite so interesting).
%
% The parameters were guessed by trial and error.
%
% Parameters:
%
%   filename     The name to give to the animation avi file.
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   compress     Boolean specifying whether to use lossy compression to
%                  (significantly) reduce the file size.  The degree of
%                  compression can be modified by changing the source code.
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 7/15/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Some animation parameters.
frames = 101;
deleteLabels = 0;               % Labels tend to float around and look bad.
qualityValue = 90;		% Even 100 gets significant compression.

%---------------------------------------------------------------------------
% Speed of motion normal to the interface.
aValue = 0.20;

% Speed of rotation (radians per unit time).
rotation = -0.75 * pi;

% Rigid body rotation or not?
rigid = 0;

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.0;                  % End time.
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (frames - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Visualize the 3D reachable set.
displayType = 'contour';

%---------------------------------------------------------------------------
% Use periodic boundary conditions?
periodic = 0;

% Create the grid.
g.dim = 2;
g.min = -1;
g.dx = 1 / 100;
if(periodic)
  g.max = (1 - g.dx);
  g.bdry = @addGhostPeriodic;
else
  g.max = +1;
  g.bdry = @addGhostExtrapolate;
end
g = processGrid(g);

%---------------------------------------------------------------------------
% What kind of display?
if(nargin < 3)
  switch(g.dim)
   case 1
    displayType = 'plot';
   case 2
    displayType = 'contour';    
   case 3
    displayType = 'surface';
   otherwise
    error('Default display type undefined for dimension %d', g.dim);
  end
end

%---------------------------------------------------------------------------
% Create initial conditions (star shaped interface centered at origin).
%   Note that in the periodic BC case, these initial conditions will not be
%   continuous across the boundary.  Regardless of boundary conditions, this
%   initial function will be far from signed distance (although it is
%   definitely an implicit surface function).  In practice, we'll just
%   ignore these little details.
points = 7;
shift = 2.5;
scale = 0.20;
data = zeros(size(g.xs{1}));
[ theta, r ] = cart2pol(g.xs{1}, g.xs{2});
data = r - scale * (cos(points * theta) + shift);
data0 = data;

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'medium';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
%   Same accuracy is used by both components of motion.
switch(accuracy)
 case 'low'
  derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Set up motion in the normal direction.
normalFunc = @termNormal;
normalData.grid = g;
normalData.speed = aValue;
normalData.derivFunc = derivFunc;

%---------------------------------------------------------------------------
% Default rotation attempts to recreate O&F figure.
if(nargin < 2)
  rigid = 0;
end

% Create linear flow field representing rigid rotation.
linearA = rotation * [ 0 1; -1 0 ];
linearV = cellMatrixMultiply(num2cell(linearA), g.xs);

rotationFunc = @termConvection;
rotationData.grid = g;
rotationData.derivFunc = derivFunc;

if(rigid)
  rotationData.velocity = linearV;
else
  % If rotation is not rigid, slow it down at smaller radii.
  % We already have polar coordinates from initial condition calculation
  rotationData.velocity = cellMatrixMultiply(r.^2, linearV);
end

%---------------------------------------------------------------------------
% Combine components of motion.
schemeFunc = @termSum;
schemeData.innerFunc = { normalFunc; rotationFunc };
schemeData.innerData = { normalData; rotationData };

%---------------------------------------------------------------------------
% Initialize Display
fig = figure;
h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);
axis(g.axis);

%---------------------------------------------------------------------------
% Now establish the animation stuff.

% Set up the figure window nicely (you can change the resolution here).
set(fig, 'Position', [ 100 100 384 384 ], 'color', 'white');
box on

% Turn off stuff we don't want -- animations are typically too busy looking.
if(deleteLabels)
  set(gca, 'XTick', [], 'YTick', []);
  delete(get(gca, 'Title'));
  delete(get(gca, 'XLabel'));
  delete(get(gca, 'YLabel'));
end

% Create the avi file (choose a smaller qualityValue to get smaller files).
if((nargin < 3) | compress)
  mov = avifile(filename, 'quality', qualityValue);
else
  mov = avifile(filename, 'compression', 'none');
end

%---------------------------------------------------------------------------
% Capture the first frame.

frame = 1;
drawnow;
f = getframe(gcf);
mov = addframe(mov, f);

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  % Delete last visualization if necessary.
  delete(h);

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
  axis(g.axis);

  % Get next frame.
  frame = frame + 1;
  drawnow;
  f = getframe(gcf);
  mov = addframe(mov, f);

end

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);

% We're finished with the movie.
mov = close(mov);
