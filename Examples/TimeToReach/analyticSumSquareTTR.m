function mttr = analyticSumSquareTTR(radius, grid)
% analyticSumSquareTTR: analytic solution special holonomic time to reach.
%
%   mttr = analyticSumSquareTTR(radius, grid)
%  
% Computes the analytic minimum time to reach each node in the grid
%   for the holonomic 2D integrator under unit bounded input.
%
% This routine is specialized to compute the analytic solution for
%   a square target set when the input is bounded in 1 norm.
%
% The dynamics are
%
%    \dot x    = b_1
%	 \dot y    = b_2
%
%   where input ||b||_1 \leq 1 is trying to hit the target.
%
% The true solution's level sets look like octagons with long diagonal faces.
%
% Parameters:
%
%   radius       Half length of the sides of the target square.
%                  Target square is centered at the origin.
%   grid         Grid structure on which data is to be computed.
%
%   mttr         Minimum time to reach function.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 12/07/04

  % Time to reach the origin when input is 1-norm bounded.
  mttr = abs(grid.xs{1}) + abs(grid.xs{2});

  % In fact, from most initial states we only have to reach the corners
  %   of the target square.
  mttr = mttr - 2 * radius;

  % For states directly in line with the initial square, it is even easier.
  mttr = chopCorner(mttr, grid.xs{1}, grid.xs{2}, radius);
  mttr = chopCorner(mttr, grid.xs{2}, grid.xs{1}, radius);

  % Finally, make sure we do not have negative time to reach.
  mttr = max(0, mttr);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function new = chopCorner(old, constrGrid, distGrid, radius)
% chopCorner: chop off a corner of a diamond.
%
%   new = chopCorner(old, constrGrid, distGrid, radius)
%
% This function takes a diamond shaped time to reach function, and
%   chops off the corners to make it octagonal.
%
% Parameters:
%   old          Initial time to reach function.
%   constrGrid   State node array for the dimension which constrains
%                  the width of the chopping.
%   distGrid     State node array for the dimension which produces
%                  the replacement data.
%   radius       Controls the width and offset of the replacement data.
%
%   new          Resultant time to reach function.
%
% Ian Mitchell 12/07/04

  % Copy the time to reach function.
  new = old;

  % Find the nodes which should be modified.
  inSquare = find(abs(constrGrid) < radius);

  % Modify those nodes.
  new(inSquare) = abs(distGrid(inSquare)) - radius;
