function [ data, g, data0 ] = ...
           reinitCircle(apply_fix, accuracy, show_nodes, grid_anisotropy, dim) 
% reinitCircle: demonstrate the subcell reinitialization fix on a circle
%
% [ data, g, data0 ] = ...
%          reinitCircle(apply_fix, accuracy, show_nodes, grid_anisotropy, dim)
%
% Demonstrates the fix specified in
%
%   Giovanni Russo & Peter Smereka, "A Remark on Computing Distance
%   Functions," J. Computational Physics, v. 163, pp. 51-67 (2000),
%   doi:10.1006/jeph.2000.6553
%
% for nodes near the interface of the reinitialization equation.  The
% example is taken from section 4, and we attempt here to recreate figures 7
% and 8.  Note that the subcell fix is applied only to nodes near the
% interface.  The rest of the reinitialization machinery is the standard
% ToolboxLS methods (see termReinit for details), so the results may differ
% from those in Russo & Smereka.
%
% In test runs, it appears that the recreation of figure 7 shows less but
% still significant volume loss when compared to the original figure 7.  The
% recreation of figure 8 is essentially identical to the original figure 8.
%
% Note that visualization is done here using Matlab's standard contour command,
% while Russo & Smereka do not specify how they generate their plots.
%
% This function was originally designed as a script file, so most of the
% options can only be modified in the file.  For example, edit the file to
% change the grid dimension, boundary conditions, flow field parameters,
% etc.
%
% Input Parameters:
%
%   apply_fix: Boolean.  Should Russo & Smereka's fix be applied?
%   Optional.  Default = True.
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
%   show_nodes: Boolean.  Plot the location of the grid nodes as well?
%   The grid is coarse (17 x 17), and plotting the nodes shows clearly
%   that reinitialization without the subcell fix allows the interface to
%   move across several nodes.  Optional.  Default = False.
%
%   grid_anisotropy: Double.  Ratio of horizontal dx to vertical dx.  The
%   vertical dx is always fixed at 10/16 as in R&S (in higher dimensions,
%   all other dx are the same as vertical dx).  A ratio of 2 will result in
%   only 8 horizontal grid cells instead of the normal 16.  Can be used to
%   demonstrate that the subcell fix works on grids whose cells are not
%   square.  Values less than one are permitted, but will result in a
%   smaller CFL timestep restriction and hence the iteration counts will not
%   be directly comparable to those in R&S.  Optional.  Default = 1.
%
%   dim: Integer in [ 2, 3 ].  Dimension in which to run the simulation.
%   Values greater than 3 should work, but there is no way to visualize the
%   results.  Optional.  Default = 2.
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
% Ian Mitchell, 5/16/07

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
if(nargin < 1)
  apply_fix = true;
end

if(nargin < 2)
  accuracy = 'low';
end

if(nargin < 3)
  show_nodes = false;
end

if(nargin < 4)
  grid_anisotropy = 1;
end

if(nargin < 5)
  dim = 2;
end

%---------------------------------------------------------------------------
% Integration parameters.

% Which iterations to plot
plot_iter = [ 0; 160; 320; 480; 640; 800 ];

% Choose a large enough time interval that we will achieve all of the
% iterations we want to plot.  We will quit integrating when we reach the
% maximum iteration, so tMax should be chosen huge.
t0 = 0;
tMax = max(plot_iter) * 100;

% Pause after each plot?
if(dim ~= 2)
  pauseAfterPlot = 1;
else
  pauseAfterPlot = 0;
end

% Some 3D visualization settings.
facecolor_next = 'red';
facecolor_last = 'blue';
alpha_last = 0.4;
alpha_next = 1.0;
if(dim == 3)
  disp('Solid red is current timestep; transparent blue is last timestep');
end

%---------------------------------------------------------------------------
% Create the grid.
g.dim = dim;
g.min = -5;
g.max = +5;
g.dx = 10/16 * [ grid_anisotropy; ones(dim - 1, 1) ];
g.bdry = @addGhostExtrapolate;
g = processGrid(g);

%---------------------------------------------------------------------------
% Initial conditions -- same as (18), but in any dimension.
center = zeros(dim, 1);
radius = 4;
data = shapeSphere(g, center, radius);
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
% appropriate iterations.  For some reason, R&S use a CFL restriction of
% 0.5 (specified in the captions of figures 7 & 8), so we will too.
integratorOptions = odeCFLset('factorCFL', 0.5, 'singleStep', 'on');

%---------------------------------------------------------------------------
% Initialize Display
figure;
axis equal;
axis(g.axis);
hold on;

if show_nodes
  % Plot the nodes in a very light grey.
  h_nodes = plot(g.xs{1}, g.xs{2}, '.');
  set(h_nodes, 'MarkerEdgeColor', 0.8 * ones(3,1));
end

if(dim == 2)
  if apply_fix
    title('Approximate recreation of figure 8 of Russo & Smereka');
  else
    title('Approximate recreation of figure 7 of Russo & Smereka');
  end
end

if(any(plot_iter == 0))
  switch(dim)
   case 2
    contour(g.xs{1}, g.xs{2}, data, [ 0 0 ], 'b-');
   case 3
    h_delete = [];
    h_last = patch(isosurface(g.xs{:}, data, 0));
    set(h_last, 'FaceColor', facecolor_next, 'EdgeColor', 'none', ...
                'FaceAlpha', alpha_next);
    lighting phong;
    camlight right;
    view(3)
   otherwise
    warning('Unable to visualize results in dimensions greater than 3.'); %#ok<WNTAG>
  end
  drawnow;
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

  % If we want to see this iteration.
  if(any(iteration == plot_iter))

    if(pauseAfterPlot)
      % Wait for last plot to be digested.
      pause;
    end

    switch(dim)
     case 2
      contour(g.xs{1}, g.xs{2}, data, [ 0 0 ], 'b-');
     case 3
      delete(h_delete);
      h = patch(isosurface(g.xs{:}, data, 0));
      set(h_last, 'FaceColor', facecolor_last, 'EdgeColor', 'none', ...
                'FaceAlpha', alpha_last);
      set(h, 'FaceColor', facecolor_next, 'EdgeColor', 'none', ...
             'FaceAlpha', alpha_next);
      h_delete = h_last;
      h_last = h;
      lighting phong;
      view(3)
    end
    drawnow
  end  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);
