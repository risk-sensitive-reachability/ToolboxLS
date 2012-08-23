% ellipseError: script to recreate figures 10 & 11 from Russo & Smereka.
%
% This script demonstrates that the subcell fix from
%
%   Giovanni Russo & Peter Smereka, "A Remark on Computing Distance
%   Functions," J. Computational Physics, v. 163, pp. 51-67 (2000),
%   doi:10.1006/jeph.2000.6553
%
% is convergent as the grid is refined, and that reaches an equilibrium
% where the interface (zero level set) stops moving.  The script recreates
% figures 10 & 11 from the R&S paper, by repeated runs of reinitEllipse.m
%
% Note that the subcell fix is applied only to nodes near the interface.
% The rest of the reinitialization machinery is the standard ToolboxLS
% methods (see termReinit for details), so the results may differ from those
% in Russo & Smereka.
%
% In experimental runs, it appears that the standard ToolboxLS
% reinitialization procedure is better than that used in R&S, in the sense
% that the rate of divergence of the error with respect to time is lower.
% The errors appearing in the recreated figure 11 are more than an order of
% magnitude off (although they show the same general behaviour), leading me
% to believe that either (27) is incorrect or that I have implemented it
% incorrectly.

% Copyright 2007 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 5/16/07

accuracy = 'low';
%accuracy = 'medium';

subcell_grids = [ 50, 100, 200 ];
no_subcell_grids = [ 200 ];

subcell_data = cell(length(subcell_grids), 1);
no_subcell_data = cell(length(no_subcell_grids), 1);

f10 = figure;
hold on;
f11 = figure;
hold on;

for i = 1 : length(subcell_grids)
  fprintf('\nReinit with subcell fix grid size %d: ', subcell_grids(i));
  subcell_data{i} = reinitEllipse(1, accuracy, subcell_grids(i), 0);
  figure(f10);
  h = plot(subcell_data{i}.t, subcell_data{i}.distance_error);
  set(h, 'LineWidth', 2, 'LineStyle', '-');
  figure(f11);
  h = plot(subcell_data{i}.t, subcell_data{i}.front_error);
  set(h, 'LineWidth', 2, 'LineStyle', '-');
end

for i = 1 : length(no_subcell_grids)
  fprintf('\nReinit without subcell fix grid size %d: ', no_subcell_grids(i));
  no_subcell_data{i} = reinitEllipse(0, accuracy, no_subcell_grids(i), 0);
  figure(f10);
  h = plot(no_subcell_data{i}.t, no_subcell_data{i}.distance_error);
  set(h, 'LineWidth', 2, 'LineStyle', '--');
  figure(f11);
  h = plot(no_subcell_data{i}.t, no_subcell_data{i}.front_error);
  set(h, 'LineWidth', 2, 'LineStyle', '--');
end


figure(f10);
set(gca, 'YScale', 'log');
axis([ 0, 15, 1e-2, 1e4]);
figure(f11);
set(gca, 'YScale', 'log');
axis([ 0, 15, 1e-5, 2e0 ]);
