function [ data, g, data0 ] = animateDumbbell(accuracy)
% animateDumbbell: animate figure 14.2 from Sethian
%
%   [ data, g, data0 ] = animateDumbbell(accuracy)
%
% Essentially the same as dumbbell1, but creates avi file.
%
% Animates the shrinking dumbbell from figure 14.2 in Sethian
%   showing motion by mean curvature of a 3D dumbbell shaped region.
%   This example is interesting because it shows pinch off and separation
%   of the implicit surface.
%
% The grid parameters have been chosen to try and get close to the figures.
%   The user can specify the curvature multiplier within the file.
%  
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters:
%
%   accuracy     Controls the order of the time approximation.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%                Note that this parameter has no effect on the order
%                  of the spatial approximation, which is always order 2.
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 5/17/04

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
deleteLabels = 1;               % Labels tend to float around and look bad.
qualityValue = 99;		% Even 100 gets significant compression.
filename = 'dumbbell';

%---------------------------------------------------------------------------
% Default values.
if(nargin < 1)
  accuracy = 'medium';
end

% Curvature multiplier.
%   You can fiddle with either this or tMax to achieve the same effect.
b = 0.025;

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1;                    % End time.
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (frames - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Pause after each plot?
pauseAfterPlot = 0;

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 3;
g.min = [ -1.0; -0.5; -0.5 ];
g.max = [ +1.0; +0.5; +0.5 ];
g.dx = 1 / 50;
g.bdry = @addGhostExtrapolate;
g = processGrid(g);

%---------------------------------------------------------------------------
% Create initial conditions (a dumbbell)
radius = 0.3;			% Radius of the dumbbell spheres.
offset = 0.5;			% Offset of the center of the dumbbell spheres.
width = 0.2;			% Width of the dumbbell center cylinder.

% Right sphere.
right = sqrt((g.xs{1} - offset).^2 + g.xs{2}.^2 + g.xs{3}.^2) - radius;

% Left sphere.
left = sqrt((g.xs{1} + offset).^2 + g.xs{2}.^2 + g.xs{3}.^2) - radius;

% Center cylinder, runs horizontally to the middle of the spheres.
center = max(abs(g.xs{1}) - offset, sqrt(g.xs{2}.^2 + g.xs{3}.^2) - width);

% Union the three portions together.
data = min(center, min(left, right));
data0 = data;

%---------------------------------------------------------------------------
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
% Set up curvature motion.
schemeFunc = @termCurvature;
schemeData.grid = g;
schemeData.curvatureFunc = @curvatureSecond;
schemeData.b = b;

%---------------------------------------------------------------------------
% Set up first frame.
  axis3D = [ -0.8, +0.8, -0.3, +0.3, -0.3, +0.3 ];
  fig = figure;
  h = visualizeLevelSet(g, data, 'surface', level);
  axis equal; axis(axis3D);
  camlight right;

%---------------------------------------------------------------------------
% Now establish the animation stuff.

% Set up the figure window nicely (you can change the resolution here).
set(fig, 'Position', [ 100 100 400 300 ], 'color', 'white');

% Turn off stuff we don't want -- animations are typically too busy looking.
if(deleteLabels)
  set(gca, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', []);
  delete(get(gca, 'Title'));
  delete(get(gca, 'XLabel'));
  delete(get(gca, 'YLabel'));
  delete(get(gca, 'ZLabel'));
end

% Create the avi file (choose a smaller qualityValue to get smaller files).
mov = avifile(filename, 'quality', qualityValue);

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

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % Create new visualization.
  delete(h);
  h = visualizeLevelSet(g, data, 'surface', level);

  frame = frame + 1;
  drawnow;
  f = getframe(gcf);
  mov = addframe(mov, f);

end

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);

% We're finished with the movie.
mov = close(mov);
