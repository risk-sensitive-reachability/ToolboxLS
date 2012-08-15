function [ mttr, data, gridOut, data0 ] = convectionTTR(front, accuracy,gridIn)
% convectionTTR: demonstrate time to reach on a trivial convective flow field.
%
%   [ mttr, data, gridOut, data0 ] = convectionTTR(front, accuracy, gridIn)
%  
% In this example we calculate the minimum time to reach a front
%   under simple convective motion.  Several different options of front shapes
%   are provided.
%
% The dynamics are
%
%        \dot x    = 1
%	 \dot y    = 0
%
% The goal is slightly more ambitious than the standard reach set, since
%   we would also like to record the minimum time to reach function by
%   tracking the time at which the reach set arrives at each node.
%
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   aircraft parameters, etc.
%
% Parameters:
%
%   front        Shape of front.
%                  'sine'        Use a sine function (default).
%                  'flat'        Use a hyperplane.
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   gridIn       Computational grid.
%                  Default is dim = 2; min = -1; max = +1; N = 101;
%                  and bdry = @addGhostExtrapolate.
%
%   mttr         Minimum time to reach function at t_max.
%   data         Implicit surface function at t_max for the reach set.
%   gridOut      Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.
%                  For data0 >= 0, this is also the analytic time to reach.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 12/06/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.

%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.6;                  % End time.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% Problem Parameters.

% Vertical speed.
speed = 1;

% Should we reinitialize after every step?
reinitialize = 1;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Visualize the reachable set.
displayType = 'contour';

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;

%---------------------------------------------------------------------------
% Create the grid (if necessary).

if(nargin < 3)
  g.dim = 2;
  g.dx = 0.02;
  g.min = [ -1; -1 ];
  g.max = [ +1; +1 - g.dx ];
  g.bdry = { @addGhostExtrapolate; @addGhostPeriodic };
else
  g = gridIn;
end

g = processGrid(g);

if(nargout > 2)
  gridOut = g;
end

%---------------------------------------------------------------------------
% Create initial conditions.
if(nargin < 1)
  front = 'sine';
end

switch(front)
 case 'sine'
  data = g.xs{1} - 0.2 * sin(2 * pi * g.xs{2}) + 0.7;
 case 'flat'
  data = g.xs{1} + 0.7;
 otherwise
  error('Unknown front shape: %s', front);
end

if(nargout > 3)
  data0 = data;
end

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = { speed; 0 };
schemeData.grid = g;

%---------------------------------------------------------------------------
if(nargin < 2)
  accuracy = 'medium';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.75, 'stats', 'on');

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
% Restrict the Hamiltonian so that reachable set only grows.
%   The Lax-Friedrichs approximation scheme MUST already be completely set up.
innerFunc = schemeFunc;
innerData = schemeData;
clear schemeFunc schemeData;

% Wrap the true Hamiltonian inside the term approximation restriction routine.
schemeFunc = @termRestrictUpdate;
schemeData.innerFunc = innerFunc;
schemeData.innerData = innerData;
schemeData.positive = 0;

%---------------------------------------------------------------------------
% Set up minimum time to reach recording using postTimestepFunc.
integratorOptions = odeCFLset(integratorOptions, ...
                              'postTimestep', @postTimestepTTR);

% Initialize the minimum time to reach function by calling
%   the postTimestepFunc once with initial data.
y = data(:);
[ y, schemeData ] = feval(@postTimestepTTR, t0, y, schemeData);
data = reshape(y, g.shape);

%---------------------------------------------------------------------------
if(reinitialize)

  % We need some more parameters to drive the reinitialization.
  %   But most can just be the defaults.
  schemeData.grid = g;
  schemeData.reinitAccuracy = accuracy;

  % Register the entire cell vector of postTimestep routines.
  integratorOptions = odeCFLset(integratorOptions, 'postTimestep', ...
                                { @postTimestepReinit; @postTimestepTTR });
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
  [ t, y, schemeData ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
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
fprintf('Total execution time %g seconds\n', endTime - startTime);

% Extract the minimum time to reach function from the schemeData structure.
%   Reshape it into an array, and replace the +inf entries with NaN to
%   make the visualization more pleasant.
mttr = reshape(schemeData.ttr, g.shape);
unreached = find(isinf(mttr));
mttr(unreached) = NaN;
