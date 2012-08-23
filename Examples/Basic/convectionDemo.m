function [ data, g, data0 ] = convectionDemo(flowType, accuracy, displayType)
% convectionDemo: demonstrate a simple convective flow field.
%
%   [ data, g, data0 ] = convectionDemo(flowType, accuracy, displayType)
%  
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters:
%
%   flowType     String to specify type of flow field.
%                  'constant'    Constant flow field xdot = k (default).
%                  'linear'      Linear flow field xdot = A x.
%                  'constantRev' Constant flow field, negate at t_half.
%                  'linearRev'   Linear flow field, negate at t_half.
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst (default).
%                  'medium'      Use odeCFL2 and upwindFirstENO2.
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
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
% Ian Mitchell, 2/9/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

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
% Most of the time in constant flow case, we want flow in a 
%   distinguished direction, so assign first dimension's flow separately.
constantV = 0 * ones(g.dim);
constantV(1) = 2;
constantV = num2cell(constantV);

% Create linear flow field xdot = A * x
linearA = 2 * pi * [ 0 1 0 0; -1 0 0 0; 0 0 0 0; 0 0 0 0 ];
%linearA = eye(4);
indices = { 1:g.dim; 1:g.dim };
linearV = cellMatrixMultiply(num2cell(linearA(indices{:})), g.xs);

%---------------------------------------------------------------------------
if(nargin < 1)
  flowType = 'constant';
end

% Choose the flow field.
switch(flowType)
  
 case 'constant'
  v = constantV;
  
 case 'linear'
  v = linearV;
  
 case 'constantRev'
  v = @switchValue;
  schemeData.one = constantV;
  schemeData.two = cellMatrixMultiply(-1, constantV);
  schemeData.tSwitch = 0.5 * tMax;
 
 case 'linearRev'
  v = @switchValue;
  schemeData.one = linearV;
  schemeData.two = cellMatrixMultiply(-1, linearV);
  schemeData.tSwitch = 0.5 * tMax;

 otherwise
  error('Unknown flowType %s', flowType);
  
end

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
% Create initial conditions (a circle/sphere).
%   Note that in the periodic BC case, these initial conditions will not
%   be continuous across the boundary unless the circle is perfectly centered.
%   In practice, we'll just ignore that little detail.
center = [ -0.4; 0.0; 0.0; 0.0 ];
radius = 0.35;
data = zeros(size(g.xs{1}));
for i = 1 : g.dim
  data = data + (g.xs{i} - center(i)).^2;
end
data = sqrt(data) - radius;
data0 = data;

%---------------------------------------------------------------------------
if(nargin < 2)
  accuracy = 'low';
end

% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = v;
schemeData.grid = g;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  schemeData.derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  schemeData.derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  schemeData.derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  schemeData.derivFunc = @upwindFirstWENO5;
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
  [ figure_az, figure_el ] = view;

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

  % Restore view.
  view(figure_az, figure_el);
  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function out = switchValue(t, data, schemeData) %#ok<INUSL>
% switchValue: switches between two values.
%
%  out = switchValue(t, data, schemeData)
%
% Returns a constant value:
%           one     for t <= tSwitch;
%           two     for t >  tSwitch.
%
% By setting one and two correctly, this function can implement
%   the velocityFunc prototype for termConvection;
%   the scalarGridFunc prototype for termNormal, termCurvature and others;
%   and possibly some other prototypes...
%
% Parameters:
%   t            Current time.
%   data         Level set function.
%   schemeData   Structure (see below).
%
%   out          Either schemeData.one or schemeData.two.
%
% schemeData is a structure containing data specific to this type of 
%   term approximation.  For this function it contains the field(s)
%
%   .one         The value to return for t <= tSwitch.
%   .two         The value to return for t >  tSwitch.
%   .tSwitch     The time at which the switch between flow fields occurs.
%
% schemeData may contain other fields.

  checkStructureFields(schemeData, 'one', 'two', 'tSwitch');
  
  if(t <= schemeData.tSwitch)
    out = schemeData.one;
  else
    out = schemeData.two;
  end
