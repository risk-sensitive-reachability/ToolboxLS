function [ data, g, data0 ] = curvatureSpiralDemo(accuracy,initial,displayType)
% curvatureSpiralDemo: demonstrate motion by mean curvature on spiral.
%
%   [ data, g, data0 ] = curvatureSpiralDemo(accuracy, initial, displayType)
%
% Recreates figure 4.1 from O&F chapter 4, showing motion by mean curvature
%   of a spiral interface in two dimensions.  The spiral is constructed 
%   in one of two ways:
%
%   a) From an ellipse in polar coordinates, whose parameters were 
%      guessed by trial and error.
%   b) From a collection of points lying on the ellipse (requires
%      iterative generation of implicit surface function).
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
%                Note that the spatial approximation is always second order.
%                  'low'         Use odeCFL1.
%                  'medium'      Use odeCFL2 (default).
%                  'high'        Use odeCFL3.
%   initial      Controls how to generate the initial implicit surface fcn.
%                  'ellipse'     From an ellipse (default).
%                  'points'      From points on the surface.
%   displayType  String to specify how to display results.
%                  The specific string depends on the grid dimension;
%                  look at the subfunction visualize to see the options
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
% Ian Mitchell, 2/13/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Curvature speed parameter.
b = 0.075;

% How should the initial data be generated?
if(nargin < 2)
  initialDataFrom = 'ellipse';
else
  initialDataFrom = initial;
end

% Reinitialize to get signed distance function at the beginning?
reinitToStart = 0;

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.0;                  % End time.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 1000 * eps;

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
% Create initial conditions (spiral around the origin).
%   Note that in the periodic BC case, these initial conditions will not be
%   continuous across the boundary.  

switch(initialDataFrom)
 case 'ellipse'
  % Regardless of boundary conditions, this initial function will be far from
  %   signed distance (although it is definitely an implicit surface
  %   function). If you want a well behaved signed distance function,
  %   reinitialize.

  % Basic ellipsoid's parameters
  M = diag([ 0.06, 2 * pi ].^(-2));
  d = 1;
  offsetR = 0.5;
  offsetTh = -pi/3;
  
  % To get a spiral, rotate the ellipsoid slightly.
  angleRot = pi / 100;
  matrixRot = [ cos(angleRot), -sin(angleRot); sin(angleRot), cos(angleRot) ];
  angleM = matrixRot * M * matrixRot';
  
  % Now use the home built spiral generator.
  data = spiralFromEllipse(g, angleM, d, offsetR, offsetTh);

 case 'points'
  samples = 200;
  data = spiralFromPoints(g, samples);

 otherwise
  error('Unknown initial string %s', initialDataFrom);
end

data0 = data;

% Reinitialize if requested (with same level of accuracy as main computation).
if(reinitToStart)
  % The maximum travel of the reinitialization wavefront should only be
  %   about a quarter of the grid size.
  tMaxReinit = 0.25 * norm(g.max - g.min);
  % We're willing to quit early if the results look good.
  errorMax = 0.01;
  data = signedDistanceIterative(g, data, tMaxReinit, errorMax, accuracy);
end

%---------------------------------------------------------------------------
% Set up motion by mean curvature with constant coefficient.
schemeFunc = @termCurvature;
schemeData.grid = g;
schemeData.curvatureFunc = @curvatureSecond;
schemeData.b = b;

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'medium';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  integratorFunc = @odeCFL1;
 case 'medium'
  integratorFunc = @odeCFL2;
 case 'high'
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

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
