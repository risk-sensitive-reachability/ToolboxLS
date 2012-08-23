% exerciseKP529: Solve Exercise 17.1.2 from Kloeden & Platen, pp.529-531
%
% This script solves Exercise 17.1.2, pp.529-531 from 
%   Kloeden & Platen, "Numerical Solution of SDEs", 3rd edition.
%
% The terminal value PDE is
%
%       D_t u + 0.5 * laplacian(u) + V(x) * u = 0
%
%   with V(x) = 0.5 * x^2 and terminal conditions u(T,x) = 1.
%
% We will solve forward in time by negating the two spatial terms.
%
% The analytic solution quoted there is (after substituting t := T - t):
%
%       u(t,x) = exp[ 0.5 * t
%                     + 0.5 * x^2 * (1 - exp(2*t))/(1 + exp(2*t))
%                     + 0.5 * ln(2 / (1 + exp(2*t))) ]
%
% Problem parameters can be modified in the file.  For example, edit the 
%   file to change the grid dimension, boundary conditions, etc.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 8/26/04

%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 0.5;                  % Final time.
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
%   Use linear extrapolation toward zero (for lack of anything better).
g.dim = 1;
g.min = -2;
g.dx = 1 / 50;
g.max = +2;
%g.bdry = @addGhostExtrapolate;
g.bdry = @addGhostExtrapolate2;
g.bdryData.towardZero = 1;
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
% What is our dependent variable, 
%   and what grid nodes should we examine for results?
switch(g.dim)
 case 1
  % Use the only dimension we are given.
  indexX = 1 : g.N(1);
  stateX = g.xs{1};
 case 2
  % The operative dimension is the diagonal x_1 = x_2.
  indexX = 1 : g.N(1) + 1 : prod(g.N);
  stateX = (g.xs{1} + g.xs{2}) * invSqrt2;
 otherwise
  error('Dimension %d not currently allowed.', dim);
end

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
% Create initial conditions.
data = ones(g.shape);
data0 = data;

%---------------------------------------------------------------------------
accuracy = 'medium';

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
% Set up scheme for the discounting term.
%   Need to negate this term, since we are solving a terminal value PDE.
%   Note that original statement in KP is incorrect; V(t,x) = -0.5 * x^2.
discountFunc = @termDiscount;
discountData.grid = g;
discountData.lambda = 0.5 * stateX.^2;

% Set up scheme for the second derivative term
%   This method computes a Laplacian, although it is not terribly efficient.
%   No need to negate, since this term's toolbox prototype is already negated.
laplacianFunc = @termTraceHessian;
laplacianData.grid = g;
laplacianData.L = 0.5 * eye(g.dim);
laplacianData.R = eye(g.dim);
laplacianData.hessianFunc = @hessianSecond;

% Combine the spatial approximation schemes.
schemeFunc = @termSum;
schemeData.innerFunc = { discountFunc; laplacianFunc };
schemeData.innerData = { discountData; laplacianData };

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

%---------------------------------------------------------------------------
% Compare with analytic solution.

% State along the line we are examining.
lineX = stateX(indexX);

% Where to put the legend (according to help entry for legend)?
legendLocation = 0;

% Analytic solution.
analytic = exp(0.5 * tMax ...
               + 0.5 * lineX.^2 * (1 - exp(2*tMax)) / (1 + exp(2*tMax)) ...
               + 0.5 * log(2 / (1 + exp(2*tMax))));

% Generate a comparison figure.
figure;
subplot(2,1,1);
plot(lineX, data(indexX), 'b-o', lineX, analytic, 'r-x');
legend('from PDE', 'analytic', legendLocation);
xlabel('x');  ylabel('u(x)');
title([ 'Solution to PDE at t = ' num2str(tMax) ]);

subplot(2,1,2);
semilogy(lineX, abs(analytic - data(indexX)), 'k-*');
xlabel('x');  ylabel('error');
