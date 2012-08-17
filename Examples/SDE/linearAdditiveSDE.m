function [ data, g, data0 ] = linearAdditiveSDE(payoff, a, b, tf, dim,accuracy)
% linearAdditive: demonstrate linear SDE with additive noise
%
%   [ data, g, data0 ] = linearAdditiveSDE(payoff, a, b, tf, dim, accuracy)
%
% This file computes the expected payoff for motion according to
%   the linear ODE with additive noise
%
%	dx = a * x * dt + b dW
%
%   where constants a,b are specified.
%
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters:
%
%   payoff       Expectation of what function should be computed?
%                    'x'           Compute expected terminal state (default).
%                    'x^2'         Compute expected terminal x^2.
%                                    Used to find variance.
%   a            Deterministic component of flow.  Scalar, defaults to 1.0.
%   b            Stochastic component of flow.  Scalar, defaults to 0.1.
%   tf           Final time.  Defaults to 1.0.
%   dim          Dimension of grid in which to run.  Defaults to 1.
%   accuracy     Controls the order of approximation.
%                    'low'         Use odeCFL1 and upwindFirstFirst.
%                    'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                    'high'        Use odeCFL3 and upwindFirstENO3.
%                    'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%                  Note that the Hessian term is always second order accurate.
%
%   data         Expected value of payoff at terminal time.
%   g            Grid structure on which data was computed.
%   data0        Exact payoff at terminal time.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 8/25/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
if(nargin < 4)
  tMax = 1.0;
else
  tMax = tf;
end

plotSteps = 26;              % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% Useful constant.
invSqrt2 = 1.0 / sqrt(2);

%---------------------------------------------------------------------------
% Pause after each plot?
pauseAfterPlot = 0;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;

% Delete previous plot unless we are using separate subplots.
deleteLastPlot = ~useSubplots;

%---------------------------------------------------------------------------
% Create the grid.
if(nargin < 5)
  g.dim = 1;
else
  g.dim = dim;
end
g.min = -2;
g.dx = 1 / 50;
g.max = +2;
%g.bdry = @addGhostExtrapolate;
g.bdry = @addGhostExtrapolate2;
g = processGrid(g);

%---------------------------------------------------------------------------
% We want to display the entire function, not just a level set.
  switch(g.dim)
   case 1
    displayType = 'plot';
   case 2
    displayType = 'surf';
   otherwise
    error('Display type undefined for dimension %d', g.dim);
  end

% Pass level = [] to visualizeLevelSet, so as to disable warnings.
%   See visualizeLevelSet docs for details.
level = [];

%---------------------------------------------------------------------------
% Define state x for 1D and 2D flows.
%   2D is just a 1D flow rotated 45 degrees (x = sqrt(2) * (x_1 + x_2)).

switch(g.dim)
 case 1
  stateX = g.xs{1};
 case 2
  stateX = (g.xs{1} + g.xs{2}) * invSqrt2;
 otherwise
  error('state is undefined for dimension %d', g.dim);
end

%---------------------------------------------------------------------------
% Flow parameters.

% Dynamics.
if(nargin < 2)
  a = 1;
elseif(numel(a) ~= 1)
  error('Only spatially constant flow is permitted; a must be a scalar.');
end
if(nargin < 3)
  b = 0.1;
elseif(numel(b) ~= 1)
  error('Only spatially constant flow is permitted; b must be a scalar.');
end

% What kind of dynamics are actually present?
deterministicComponent = (a ~= 0);
stochasticComponent = (b ~= 0);

switch(g.dim)
 case 1
  % Motion in the only dimension possible.
  v{1} = -a * stateX;
  L = 0.5 * b^2;
  R = 1;

 case 2
  % Motion along the diagonal x_1 + x_2.
  %   Command deal copies the same matrix into both elements of cell vector v.
  [ v{1:2} ] = deal(-a * invSqrt2 * stateX);
  L = 0.5 * (invSqrt2 * b)^2 * ones(2,2);
  R = eye(2);

 otherwise
  error('Dynamics undefined for dimension %d', g.dim);

end

%---------------------------------------------------------------------------
% Create initial conditions.
%   These encode the final payoff, whose expectation the PDE will compute.

if(nargin < 1)
  payoff = 'x';
end

switch(payoff)
 case 'x'
  data = stateX;
 case 'x^2'
  data = stateX .^ 2;
 otherwise
  error('Unknown terminal payoff choice %s', payoff);
end

data0 = data;

%---------------------------------------------------------------------------
if(nargin < 6)
  accuracy = 'medium';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  upwindDerivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  upwindDerivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  upwindDerivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  upwindDerivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Set up spatial approximation scheme for the deterministic motion.
if(deterministicComponent)
  deterministicFunc = @termConvection;
  deterministicData.grid = g;
  deterministicData.velocity = v;
  deterministicData.derivFunc = upwindDerivFunc;
end

% Set up spatial approximation scheme for the stochastic motion.
if(stochasticComponent)
  stochasticFunc = @termTraceHessian;
  stochasticData.grid = g;
  stochasticData.L = L;
  stochasticData.R = R;
  stochasticData.hessianFunc = @hessianSecond;
end

% Combine the spatial approximation schemes.
if(deterministicComponent)
  if(stochasticComponent)
    schemeFunc = @termSum;
    schemeData.innerFunc = { deterministicFunc; stochasticFunc };
    schemeData.innerData = { deterministicData; stochasticData };
  else
    schemeFunc = deterministicFunc;
    schemeData = deterministicData;
  end
else
  if(stochasticComponent)
    schemeFunc = stochasticFunc;
    schemeData = stochasticData;
  else
    error('System has zero dynamics');
  end
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
