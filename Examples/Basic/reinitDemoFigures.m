% reinitDemoFigures: figures from reinitDemo
%
% This is a script file to generate the figures for reinitDemo that
% appear in the Toolbox documentation.

% Copyright 2007 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 5/21/07

% Get the data by running reinitDemo
[ data, g, data0 ] = reinitDemo('star', 'medium', 'surf');

% The figure produced at the end shows the final level set, but we need
% to make it match a figure showing the initial conditions.  We also need
% to remove the title.
f2 = gcf;
title('');

% Create a figure with the initial conditions.
f1 = figure;
surf(g.xs{1}, g.xs{2}, data0);
axis image;
initial_axis = axis;

% Adjust the final implicit surface function to use the same axis.
figure(f2);
axis(initial_axis);

% Now a figure to show the effects of reinitialization on the zero contour.
f3 = figure;
[ ~, h_contour_final ] = contour(g.xs{1}, g.xs{2}, data, [ 0 0 ], 'b-');
hold on;
[ ~, h_contour_init ] = contour(g.xs{1}, g.xs{2}, data0, [ 0 0 ], 'r--');
axis equal;
axis(g.axis);
set(h_contour_final, 'LineWidth', 2);
set(h_contour_init, 'LineWidth', 2);
legend([ h_contour_final(1), h_contour_init(1) ], ...
       'after reinitialization', 'before reinitialization');

% And finally a figure to show that the magnitude of the gradient is
% one.  Note that Matlab's grad function does centered differences, while
% reinitialization is based on the upwind differences in the Toolbox.
% But the fact that the gradient magnitude is well behaved even for a
% different finite difference is good evidence.
f4 = figure;
[ dx0 dy0 ] = gradient(data0, g.dx(1), g.dx(2));
[ dx, dy ] = gradient(data, g.dx(1), g.dx(2));
mag0 = sqrt(dx0.^2 + dy0.^2);
mag = sqrt(dx.^2 + dy.^2);
h_mag_init = plot(sort(mag0(:)), 'r--');
hold on;
h_mag_final = plot(sort(mag(:)), 'b-');
set(h_mag_init, 'LineWidth', 2);
set(h_mag_final, 'LineWidth', 2);
mag_axis = axis;
mag_axis([ 3, 4 ]) = [ 0, 5 ];
axis(mag_axis);
legend([ h_mag_final(1), h_mag_init(1) ], ...
       'after reinitialization', 'before reinitialization');
