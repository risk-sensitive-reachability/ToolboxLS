function [ data, g, data0 ] = normalStarDemo(accuracy, reverseFlow,displayType)
% normalStarDemo: demonstrate motion by surface normal on star interface.
%
%   [ data, g, data0 ] = normalStarDemo(accuracy, reverseFlow, displayType)
%
% Recreates figure 6.1 from O&F chapter 6, showing motion by surface normal
%   of a star-shaped interface in two dimensions.  The parameters were
%   guessed by trial and error.
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
%   reverseFlow  Boolean specifies whether to reverse the motion of the
%                  flow at half time (default == 0).
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
% Ian Mitchell, 3/1/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Speed of motion normal to the interface.
aValue = 0.25;

% Use the time dependent motion to reverse the flow?
if(nargin < 2)
  useTimeDependent = 0;
else
  useTimeDependent = reverseFlow;
end

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
% Set up motion in the normal direction (derivative choice is set below).
schemeFunc = @termNormal;
schemeData.grid = g;

if(useTimeDependent)
  % Time dependent flow field, switches direction.
  %   For t <= tHalf, grow outward.
  %   For t >  tHalf, shrink inward.
  schemeData.speed = @switchValue;
  schemeData.tSwitch = 0.5 * tMax;
  schemeData.one = +aValue;
  schemeData.two = -aValue;
else
  % Time independent flow field is constant.
  schemeData.speed = aValue;
end

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'medium';
end

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



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function out = switchValue(t, ~, schemeData)
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
