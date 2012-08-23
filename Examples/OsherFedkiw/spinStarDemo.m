function [ data, g, data0 ] = spinStarDemo(accuracy, rigid, displayType)
% spinStarDemo: demonstrate combined surface motion of a star shape.
%
%   [ data, g, data0 ] = spinStarDemo(accuracy, rigid, displayType)
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
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters:
%
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   rigid        Boolean specifying whether rotation should be rigid body
%                  or should recreate figure 6.2 (defaults to 0, the figure).
%   displayType  String to specify how to display results.
%                  The specific string depends on the grid dimension;
%                  look at the helper visualizeLevelSet to see the options
%                  (optional, default depends on grid dimension).
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 4/9/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Speed of motion normal to the interface.
aValue = 0.20;

% Speed of rotation (radians per unit time).
rotation = -0.75 * pi;

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.0;                  % End time.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 0;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

%---------------------------------------------------------------------------
% Use periodic boundary conditions?
periodic = 0;

% Create the grid.
g.dim = 2;
g.min = -1;
g.dx = 1 / 50;
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
f = figure;

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);

hold on;
if(g.dim > 1)
  axis(g.axis);
  daspect([ 1 1 1 ]);
end

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

  % Get correct figure, and remember its current view.
  figure(f);

  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end

  % Move to next subplot if necessary.
  if(useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);
