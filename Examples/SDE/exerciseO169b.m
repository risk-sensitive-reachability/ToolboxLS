function [ data, g, data0 ] = exerciseO169b
% exerciseO169b: Solve Exercise 8.6 from Oksendal, pp.169-170
%
% This script solves Exercise 8.6, pp.169-170 from 
%   Oksendal, "Stochastic Differential Equations", sixth edition.
%
% The initial value PDE for x \in \R is
%
%       D_t u = -\rho u + \alpha x D_x u + 0.5 \beta^2 x^2 D_x^2 u
%
%   with initial conditions u(0,x) = (x - K)^+
%   where \rho > 0, K > 0, \alpha, \beta are constant and
%   (x - K)^+ = max(x - K, 0)
%
% This PDE apparently arising in connection with deduction of the
%   Black-Scholes formula for the price of an option.
%
% The analytic solution given in the exercise is
%
%       u(t,x) = exp(-\rho / t) / \sqrt{2 \pi t}
%                \int_\R (x exp[(\alpha - 0.5 \beta^2)t + \beta y] - K)^+
%                        exp(-y^2 / 2 t) dy
%
% This analytic solution is approximated by quadrature for comparison
%   with the numerical solution for the PDE deduced by the toolkit.
%
% Based on experimentation, it appears that Matlab's quadrature routines
%   have a tough time with the max function for x values between 0 and K.
%   The problem manifests as an apparent jump in the analytic function
%   from 0 to the numerical approximation value somewhere in this x range.
%   Choosing larger quadrature widths makes it worse.
%
% Experimentation with other functions (such as x^2 instead of max(x-k,0))
%   gets the correct results, so hopefully this is a failure of the
%   quadrature rules rather than the toolbox.
%
% Problem parameters can be modified in the file.  For example, edit the 
%   file to change the grid dimension, boundary conditions, etc.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 9/9/04

%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Problem parameters
rho = 1;
alpha = 1;
beta = 1;
K = 1;

% We approximate the analytic solution with a quadrature.
%   How wide should we make the quadrature?
% If you make it too wide, experimental evidence suggests that the
%   quadrature starts to return incorrect estimates
%   (it appears to have real trouble with the max function).
quadratureWidth = 10;

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.0;                  % Final time.
plotSteps =  4;              % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% Pause after each plot?
pauseAfterPlot = 0;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

% Delete previous plot unless we are using separate subplots.
deleteLastPlot = ~useSubplots;

%---------------------------------------------------------------------------
% Create the grid.
%   Use linear extrapolation toward zero (for lack of anything better).
g.dim = 1;
g.min = -1;
g.max = +4;
g.N = 101;
%g.bdry = @addGhostExtrapolate;
g.bdry = @addGhostExtrapolate2;
g.bdryData.towardZero = 1;
g = processGrid(g);

%---------------------------------------------------------------------------
% We want to display the entire function, not just a level set.
displayType = 'plot';

% Pass level = [] to visualizeLevelSet, so as to disable warnings.
%   See visualizeLevelSet docs for details.
level = [];

%---------------------------------------------------------------------------
% Create initial conditions.
data0 = max(g.xs{1} - K, 0);
data = data0;

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
discountFunc = @termDiscount;
discountData.grid = g;
discountData.lambda = rho;

% Set up scheme for the convective term.
convectFunc = @termConvection;
convectData.grid = g;
convectData.derivFunc = upwindDerivFunc;
convectData.velocity = { -alpha * g.xs{1} };

% Set up scheme for the second derivative term.
%   This is a one dimensional problem, but the multiplier depends on state,
%   so we end up with a 1x1 cell matrix.
%   No need to negate, since this term's toolbox prototype is already negated.
laplacianFunc = @termTraceHessian;
laplacianData.grid = g;
laplacianData.L = { 0.5 * beta^2 * g.xs{1}.^2 };
laplacianData.R = 1;
laplacianData.hessianFunc = @hessianSecond;

% Combine the spatial approximation schemes.
schemeFunc = @termSum;
schemeData.innerFunc = { discountFunc; convectFunc; laplacianFunc };
schemeData.innerData = { discountData; convectData; laplacianData };

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

hL = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);
hold on;

% There is no need to plot the analytic solution at t = 0
hA = [];
axis tight;

% For documentation, make the lines thicker.
%set([ hL; hA ], 'LineWidth', 2);

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
    delete(hL);
    delete(hA);
  end

  % Move to next subplot if necessary.
  if(useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  % Create new visualization.
  hL = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow)]);
  hold on;
  hA = plotSoln(tNow, g.xs{1}, rho, alpha, beta, K, quadratureWidth);
  axis tight;

  % For documentation, make the lines thicker.
  %set([ hL; hA ], 'LineWidth', 2);
  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function h = plotSoln(t, xs, rho, alpha, beta, K, quadratureWidth)
% plotSoln: plot an approximation to the analytic solution from Oksendal
%
%   h = plotSoln(t, xs, rho, alpha, beta, K, quadratureWidth)
%
% Plots the (approximation of the) analytic solution given in Oksendal
%   by calling the analyticSoln subfunction repeatedly.
%
% Plotting is done in the current figure.
%
% Parameters:
%   t               Current time.
%   xs              Vector of points at which to plot solution.
%   rho             Problem parameter.
%   alpha           Problem parameter.
%   beta            Problem parameter.
%   K               Problem parameter.
%   quadratureWidth Width of the numerical quadrature approximation.
%
%   h               Handle to the plotted line.

solns = zeros(size(xs));
for i = 1 : length(xs);
  solns(i) = analyticSoln(t, xs(i), rho, alpha, beta, K, quadratureWidth);
end

h = plot(xs, solns, 'r--');
drawnow;


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function soln = analyticSoln(t, x, rho, alpha, beta, K, quadratureWidth)
% analyticSoln: approximate the analytic solution from Okensedal.
%
%   soln = analyticSoln(t, x, rho, alpha, beta, K, quadratureWidth)
%
% This function approximates the solution of the equation
%
%       u(t,x) = exp(-\rho / t) / \sqrt{2 \pi t}
%                \int_\R (x exp[(\alpha - 0.5 \beta^2)t + \beta y] - K)^+
%                        exp(-y^2 / 2 t) dy
%
%   from Oksendal p.170 using numerical quadrature.
%
% Parameters:
%   t               Current time.
%   x               Current state (scalar).
%   rho             Problem parameter.
%   alpha           Problem parameter.
%   beta            Problem parameter.
%   K               Problem parameter.
%   quadratureWidth Width of the numerical quadrature approximation.
%
%   soln            Approximate solution to the integral.

f = inline([ 'max(x * exp((alpha - 0.5*beta^2) * t + beta * y) - K, 0)' ...
             ' .* exp(-0.5 * y.^2 / t)' ], 'y', 't', 'x', 'alpha', 'beta','K');

integral = quad(f, -quadratureWidth, quadratureWidth, [], [], ...
                t, x, alpha, beta, K);
soln = integral * exp(-rho * t) / sqrt(2 * pi * t);
