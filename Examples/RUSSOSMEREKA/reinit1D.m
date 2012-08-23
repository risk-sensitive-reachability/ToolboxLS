function [ data, g, data0 ] = reinit1D(apply_fix, accuracy)
% reinit1D: demonstrate the subcell reinitialization fix in 1D
%
% [ data, g, data0 ] = reinit1D(apply_fix, accuracy)
%
% Demonstrates the fix specified in
%
%   Giovanni Russo & Peter Smereka, "A Remark on Computing Distance
%   Functions," J. Computational Physics, v. 163, pp. 51-67 (2000),
%   doi:10.1006/jeph.2000.6553
%
% for nodes near the interface of the reinitialization equation.  The
% example is taken from section 2 (it also appears in section 3), and we
% attempt here to recreate figures 2 and 5.  Note that the subcell fix is
% applied only to nodes near the interface.  The rest of the
% reinitialization machinery is the standard ToolboxLS methods (see
% termReinit for details), so the results may differ slightly from those in
% Russo & Smereka.
%
% Note that visualization is done here using Matlab's standard plot command,
% while the description in Russo & Smereka makes it likely that the lines in
% figures 2 and 5 are cubic spline fits to the data set.  Under the same
% restriction, if you would like to see figures 3 and 6 use the following
% command after this routine has finished running:
%
%    axis([ -0.9, +0.4, -2.5, +2.5 ]);
%
% This function was originally designed as a script file, so most of the
% options can only be modified in the file.  For example, edit the file to
% change the grid dimension, boundary conditions, flow field parameters,
% etc.
%
% Input Parameters:
%
%   apply_fix: Boolean.  Should Russo & Smereka's fix be applied?
%   Optional.  Default = 1.
%
%   accuracy: String.  Controls the order of approximations away from the
%   interface.  If the subcell fix is applied, the order of the
%   approximation near the interface is governed by the subcell fix.  The
%   choice of 'low' should correspond with Russo & Smereka's results.
%   Optional.
%
%     'low'         Use odeCFL1 and upwindFirstFirst (default).
%     'medium'      Use odeCFL2 and upwindFirstENO2.
%     'high'        Use odeCFL3 and upwindFirstENO3.
%     'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%
% Output Parameters:
%
%   data: Implicit surface function at the final iteration.
%
%   g: Grid structure on which data was computed.
%
%   data0: Implicit surface function at the beginning.

% Copyright 2007 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 5/5/07

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
%run('../addPathToKernel');

%---------------------------------------------------------------------------
if(nargin < 1)
  apply_fix = 1;
end

if(nargin < 2)
  accuracy = 'low';
end

%---------------------------------------------------------------------------
% Integration parameters.

% Which iterations to plot
plot_iter = [ 0; 3; 6; 9; 12 ];

% Choose a large enough time interval that we will achieve all of the
% iterations we want to plot.  We will quit integrating when we reach the
% maximum iteration, so tMax should be chosen huge.
t0 = 0;
tMax = max(plot_iter) * 100;

% Pause after each plot?
pauseAfterPlot = 0;

% Plotting parameters.
line_width = 2;

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 1;
g.min = -5;
g.max = +5;
g.dx = 0.5;
g.bdry = @addGhostExtrapolate;
g = processGrid(g);

%---------------------------------------------------------------------------
% Initial conditions from (9).
data = 0.5 * (g.xs{1} - 0.4 * g.dx(1)) .* (g.xs{1} + 6) + 1;
data0 = data;

%---------------------------------------------------------------------------
% Choose approximations at appropriate level of accuracy.
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

% Set up spatial approximation scheme.
schemeFunc = @termReinit;
schemeData.grid = g;
schemeData.initial = data0;
schemeData.derivFunc = derivFunc;

if apply_fix
  schemeData.subcell_fix_order = 1;
else
  schemeData.subcell_fix_order = 0;
end

% Set up time approximation scheme.  We don't bother with stats, since the
% execution is so fast.  Use single stepping so that we can plot the
% appropriate iterations.
integratorOptions = odeCFLset('factorCFL', 0.9, 'singleStep', 'on');

%---------------------------------------------------------------------------
% Initialize Display
f = figure;

% Show the zero level set and the node locations.
h = plot(g.xs{1}, zeros(size(g.xs{1})), 'r-o');
set(h, 'LineWidth', line_width);
hold on;

if apply_fix
  title('Approximate recreation of figure 5 of Russo & Smereka');
else
  title('Approximate recreation of figure 2 of Russo & Smereka');
end

if(any(plot_iter == 0))
  h = plot(g.xs{1}, data, 'b--');
  set(h, 'LineWidth', line_width);
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
iteration = 0;
while(iteration < max(plot_iter))

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, tMax ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);
  iteration = iteration + 1;
  
  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % If we want to see this iteration.
  if(any(iteration == plot_iter))
    h = plot(g.xs{1}, data, 'b-');
    set(h, 'LineWidth', line_width);
  end
  
end

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);
