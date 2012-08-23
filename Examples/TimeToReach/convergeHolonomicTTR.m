% convergeHolonomicTTR: script file to demonstrate holonomic convergence.
%
% This script file runs through a variety of grid resolutions and
%   algorithm accuracies to generate empirical convergence plots for 
%   the simple holonomic system example of minimum time to reach.
%
% Most parameters must be set within holonomicTTR, rather than here.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 12/06/04

% Solution generator (both analytic and numeric).
numeric = @holonomicTTR;

% Which norm?
whichNorm = 'sum';

% What accuracies?
accuracy = { 'low', 'medium', 'high', 'veryHigh' };
%accuracy = { 'low', 'medium' };
%accuracy = { 'high' };

% What resolutions?
%   Note that some fudging will occur (eg, add 1 to get round g.dx).
resolution = { 50; 71; 100; 141; 200; 283; 400; 566 };
%resolution = { 50; 71; 100; 141; 200 };
%resolution = { 50; 71; 100; };

err.max = zeros(length(accuracy), length(resolution));
err.two = zeros(length(accuracy), length(resolution));
err.avg = zeros(length(accuracy), length(resolution));

gridDX = zeros(length(resolution), 1);

for res = 1 : length(resolution),
  % Computational grid.
  clear g;
  g.dim = 2;
  g.min = -1.1;
  g.max = +1.1;
  g.bdry = @addGhostExtrapolate;
  
  % Increase resolution to handle increased grid size and for round g.dx.
  g.N = round(resolution{res} * 1.1 + 1);
  
  g = processGrid(g);

  gridDX(res) = min(g.dx);

  for acc = 1 : length(accuracy),

    fprintf('\nGrid Size %d, accuracy %s\n', resolution{res}, accuracy{acc});

    % Get the solutions.
    [ nttr, attr ] = feval(numeric, whichNorm, accuracy{acc}, g);

    % Compare only in interior of domain and where nttr is defined.
    interior = find((g.xs{1} >= -1) & (g.xs{1} <= +1) & ...
                    (g.xs{2} >= -1) & (g.xs{2} <= +1) & ...
                    ~isnan(nttr));

    numInterior = length(interior);

    % Collect norms.
    err.max(acc,res) = norm(nttr(interior) - attr(interior), inf);
    err.avg(acc,res) = norm(nttr(interior) - attr(interior), 1) / numInterior;
    err.two(acc,res) = norm(nttr(interior) - attr(interior), 2) ...
                          / sqrt(numInterior);

  end
end

% Plot some results.
gridSizes = cell2mat(resolution);
errorTypes = fieldnames(err);
schemeStyles = { 'm:v', 'r-.*', 'k-x', 'b--o' };
for type = 1 : length(errorTypes)
  errorType = errorTypes{type};
  f.(errorType) = figure;
  title([ 'Error analysis for ' errorType ' error.' ]);
  hold on;

  for acc = 1 : length(accuracy),
    errs = err.(errorType);
    h(acc, type) = plot(gridSizes, errs(acc,:), schemeStyles{acc});
  end

  h(length(accuracy) + 1, type) = plot(gridSizes, gridDX, 'g-');

  set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', gridSizes, ...
           'XGrid', 'on', 'XMinorGrid', 'off', ...
           'YGrid', 'on', 'YMinorGrid', 'off');
  xlabel('Grid Size');  ylabel('Error');

  % Set the axes so that they include the data (plus some space vertically).
  currAxis = axis;
  newAxis = [ min(gridSizes), max(gridSizes), ...
              10^(floor(log10(currAxis(3)))), 10^(ceil(log10(currAxis(4)))) ];
  axis(newAxis);

  % Show slopes with various orders of convergence.
  slopePoint = [ gridSizes(end-2); 100 * newAxis(3) ];
  slopeHandles = addSlopes(slopePoint, ...
                           gridSizes(end-1) - gridSizes(end-2), 'k-', ...
                           [ 0, -1, -2 ], ...
                           { '', '1', '2' });

  % Label the schemes.
  legend(h(:, type), { accuracy{:}, '\Delta x' } );

  % Make the lines and markers a bit more visible.
  set(h(:, type), 'LineWidth', 2.0, 'MarkerSize', 9.0);
end
