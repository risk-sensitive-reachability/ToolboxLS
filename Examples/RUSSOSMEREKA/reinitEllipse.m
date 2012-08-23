function [ data_out, g ] = ...
                    reinitEllipse(apply_fix, accuracy, num_nodes, do_figure_9)
% reinitCircle: demonstrate the subcell reinitialization fix on a circle
%
% [ data_out, g, data0 ] = ...
%                   reinitEllipse(apply_fix, accuracy, num_nodes, do_figure_9)
%
% Demonstrates the fix specified in
%
%   Giovanni Russo & Peter Smereka, "A Remark on Computing Distance
%   Functions," J. Computational Physics, v. 163, pp. 51-67 (2000),
%   doi:10.1006/jeph.2000.6553
%
% for nodes near the interface of the reinitialization equation.  The
% example is taken from section 4.  This file can either recreate figure 9,
% or it can be used by a driver routine to generate data for recreating
% figures 10 & 11.  Note that the subcell fix is applied only to nodes near
% the interface.  The rest of the reinitialization machinery is the standard
% ToolboxLS methods (see termReinit for details), so the results may differ
% from those in Russo & Smereka.
%
% In test runs, it appears that the recreation of figure 7 shows less but
% still significant volume loss when compared to the original figure 7.  The
% recreation of figure 8 is essentially identical to the original figure 8.
%
% Note that visualization for figure 9 is done here using Matlab's standard
% contour command, while Russo & Smereka do not specify how they generate
% their plots.
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
%   num_nodes: Integer or 2 element vector of integers.  Number of nodes in
%   the grid.  The actual number used is num_nodes + 1.  Optional.  Default
%   = 200 (the R&S choice).
%
%   do_figure_9: Boolean.  Generate a recreation of figure 9 of R&S?  If
%   not, it is assumed that this routine is being called to generate data
%   for recreating either figure 10 or 11.  Consequently, this input
%   parameter determines the length of integration and also the return
%   parameters.  Optional.  Default = 1.
%
% Output Parameters:
%
%   data_out: If do_figure_9 is true, data_out is a cell vector containing the
%   implicit surface function at each of the four iterations shown in figure
%   9.  If do_figure_9 is false, data_out is a structure with three elements,
%   each of which is a vector.
%
%     .t: Timesteps.
%
%     .distance_error: The L1 error in the signed distance function at
%     each timestep.
%
%     .front_error: The L1 error in the location of the zero level set (as
%     defined by (27)) at each timestep.
%
%   g: Grid structure on which the run was performed.

% Copyright 2007 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 5/16/07

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

% How close (relative) do we need to get to tMax to be considered finished?
small = 1000 * eps;

%---------------------------------------------------------------------------
if(nargin < 1)
  apply_fix = 1;
end

if(nargin < 2)
  accuracy = 'low';
end

if(nargin < 3)
  num_nodes = 200;
end

if(nargin < 4)
  do_figure_9 = 1;
end

%---------------------------------------------------------------------------
% Initial condition parameters from R&S (just after (24)).  Had to flip
% the sign of x0 and y0 to match figure 9.
ellipse.A = 4;
ellipse.B = 2;
ellipse.epsilon = 0.1;
ellipse.x0 = -3.5;
ellipse.y0 = -2.0;

% For some reason, R&S use a CFL restriction of 0.5 (specified in the
% caption of figure 9), so we will too.
factor_cfl = 0.5;

% Contour plotting parameters for figure 9.  The caption states that the
% contours are -1:0.2:+1.  However, looking at a cross-section through the
% contour plots at x = 0 for the lower right corner (where the solution
% should be very close to signed distance) we see five contours between +2
% and 0.  Therefore we use contours -2:0.4:+2.
contour_lines = -2:0.4:+2;

% Number of data points on the ellipse used for error evaluation (from
% the text just below figure 11).
n_sigma = 2000;

%---------------------------------------------------------------------------
% Integration parameters.

% Which iterations to plot (if we are doing figure 9).  We had to multiply
% the iteration counts shown in the caption of figure 9 by four to get
% contour plots the looked similar to those in figure 9.
plot_iter = 4 * [ 0, 10, 25, 50 ];

% This time interval should be the one used by figures 10 & 11, if we are
% not doing figure 9.
t0 = 0;
tMax = 15;

% Pause after each plot?
pauseAfterPlot = 0;

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
g.min = -5;
g.max = +5;
g.N = num_nodes + 1;
g.bdry = @addGhostExtrapolate;
g = processGrid(g);

%---------------------------------------------------------------------------
% Initial conditions from (24) with fix (y^2 in (24) should be divided by B).
ellipse.f = (ellipse.epsilon ...
             + (g.xs{1} - ellipse.x0).^2 + (g.xs{2} - ellipse.y0).^2);
data = ((sqrt((g.xs{1} / ellipse.A).^2 + (g.xs{2} / ellipse.B).^2) - 1) ...
        .* ellipse.f);
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
% appropriate iterations, or gather error data at each step.
integratorOptions = odeCFLset('factorCFL', factor_cfl, 'singleStep', 'on');

%---------------------------------------------------------------------------
% Set up the output data.
if do_figure_9
  data_out = cell(length(plot_iter), 1);
else
  % Add one step for the initial conditions.
  steps = ceil((tMax - t0) / factor_cfl) + 1;
  data_out.t = zeros(steps, 1);
  data_out.distance_error = zeros(steps, 1);
  data_out.front_error = zeros(steps, 1);
end

%---------------------------------------------------------------------------
% Initialize Display

% If we are recreating figure 9, keep track of which iteration and subplot
% we are working on.
iteration = 0;
subplot_number = 1;

if do_figure_9

  if(any(plot_iter == 0))
    % Save results.
    data_out{subplot_number} = data;
    % Add to the figure.
    f = figure;
    subplot_number = subplot_number + 1;
    % Generate the thin contours from figure 9.
    [ garbage, h_thin ] = contour(g.xs{1}, g.xs{2}, data, contour_lines, 'b-');
    hold on;
    % Generate the thick contour at zero from figure 9.
    [ garbage, h_thick ] = contour(g.xs{1}, g.xs{2}, data, [ 0 0 ], 'k-');
    set(h_thick, 'LineWidth', 2);
    title([ 'iteration ' num2str(iteration) ]);
    axis equal;
    axis(g.axis);
  end
else
  % Prepare to collect error statistics by precomputing the appropriate
  % data from (25) - (27).
  theta_p = linspace(0, 2*pi, n_sigma);
  x_p = ellipse.A * cos(theta_p);
  y_p = ellipse.B * sin(theta_p);

  % What is the distance between two neighboring points on the ellipse?  If it
  % was a circle, then this value would be constant because delta theta is
  % constant; however, on an ellipse it will vary.
  vector_x_p = [ x_p; y_p ];
  vector_x_p_offset = [ vector_x_p(:,2:end), vector_x_p(:,1) ];
  delta_x_p = sqrt(sum((vector_x_p - vector_x_p_offset).^2,1));
  
  % I suspect that vectorizing this operation would take up too much memory.
  data_true = inf * ones(g.shape);
  for p = 1 : n_sigma
    data_true = min(sqrt((g.xs{1} - x_p(p)).^2 + (g.xs{2} - y_p(p)).^2), ...
                    data_true);
  end
  data_true = data_true .* sign(data);
  
  % Collect the stats for this zeroth step.
  data_out.t(iteration + 1) = 0;
  % Error in implicit surface function (25).
  distance_error = abs(data - data_true);
  data_out.distance_error(iteration + 1) = prod(g.dx) * sum(distance_error(:));
  % Error in interface location (27).
  data_x_p = abs(interpn(g.xs{:}, data, x_p, y_p, 'cubic'));
  data_x_p_offset = [ data_x_p(2:end), data_x_p(1) ];
  front_error = sum((data_x_p + data_x_p_offset) .* delta_x_p);
  data_out.front_error(iteration + 1) = 0.5 * front_error;
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while((do_figure_9 & (iteration < max(plot_iter))) ...
      || (~do_figure_9 & (tMax - tNow > small * tMax)))

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

  if do_figure_9
    if(any(plot_iter == iteration))
      % Save results.
      data_out{subplot_number} = data;
      % Add to the figure.
      f = figure;
      subplot_number = subplot_number + 1;
      % Generate the thin contours from figure 9.
      [ garbage, h_thin ] = contour(g.xs{1}, g.xs{2}, data,contour_lines,'b-');
      hold on;
      % Generate the thick contour at zero from figure 9.
      [ garbage, h_thick ] = contour(g.xs{1}, g.xs{2}, data, [ 0 0 ], 'k-');
      set(h_thick, 'LineWidth', 2);
      title([ 'iteration ' num2str(iteration) ]);
      axis equal;
      axis(g.axis);
    end
  else
    % Collect the stats for this zeroth step.
    data_out.t(iteration + 1) = tNow;
    % Error in implicit surface function (25).
    distance_error = abs(data - data_true);
    data_out.distance_error(iteration+1) = prod(g.dx) * sum(distance_error(:));
    % Error in interface location (27).
    data_x_p = abs(interpn(g.xs{:}, data, x_p, y_p, 'cubic'));
    data_x_p_offset = [ data_x_p(2:end), data_x_p(1) ];
    front_error = sum((data_x_p + data_x_p_offset) .* delta_x_p);
    data_out.front_error(iteration + 1) = 0.5 * front_error;
  end
  
end

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);
