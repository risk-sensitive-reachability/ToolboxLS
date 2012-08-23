function [ hCurve, hAll ] = ...
                       visualizeOpenCurve(grid, curve, mask, showAll, titleStr)
% visualizeOpenCurve: plots an open curve from a vector level set.
%
%   [ hCurve, hAll ] = ...
%                      visualizeOpenCurve(grid, curve, mask, showAll, titleStr)
%
% Regular level set functions can be used only for closed curves (in 2D
% spaces).  As observed in
%
%      Peter Smereka, "Spiral Crystal Growth," 
%      Physica D 138, pp. 282-301 (2000).
%
% open curves can be represented by a vector level set:
%
%      \Gamma = \{ x | \phi(x) = 0 and \psi(x) >= 0 \}
%
% Following Smereka's method for plotting this "open curve", we define
%
%      \theta = abs(\phi), if \psi \geq 0;
%                 +1,      if \psi   <  0.
%
% and then plot the zero contour of \theta - 0.5 * max(g.dx).  This
% generates a double line that is usually visually indistinguishable from a
% single line.
%
% If the boolean parameter showAll is set, the entire zero contour of \phi
% is plotted first, and then overlapped with the \theta contour for those
% portions which are physical.
%
% Parameters:
%
%   grid         Grid structure for both the curve and mask functions.
%                  Must be two dimensional.
%   curve        Level set function whose zero level set defines the curve.
%                  In Smereka's notation, \phi.
%   mask         Level set function whose zero sublevel set defines
%                  the portion of the state space where the curve is real.
%                  In Smereka's notation, \psi.
%   showAll      Boolean specifying whether to show the entire zero contour
%                  of \phi (distinguishing between the physical step line
%                  and the artificial step line with color).  Default = 0.
%   titleStr     String to place as a title to the generated plot.
%                  Default to no title.
%
%   hCurve       Graphics handle(s) to the line(s) making up the zero 
%                  contour of the physical open curve.
%   hAll         Graphics handle(s) to the line(s) making up the zero
%                  contour of \phi.  Only available if showAll = 1; otherwise
%                  returns [].

% Copyright 2005 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/09/05

%---------------------------------------------------------------------------
% This routine only works in two dimensions.
if(grid.dim ~= 2)
  error('Open curves can only be plotted in two dimensions at present');
end

%---------------------------------------------------------------------------
% Optional arguments.
if(nargin < 4)
  showAll = 0;
end

%---------------------------------------------------------------------------
% If the user wants the entire contour.
if(showAll)
  [ garbage, hAll ] = contour(grid.xs{:}, curve, [ 0 0 ], 'r-');
  hold on;
else
  hAll = [];
end

%---------------------------------------------------------------------------
% Construct auxiliary theta function.
theta = abs(curve) .* (mask >= 0) + (mask < 0);

% Spacing between the double contour for theta.
offset = 0.5 * max(grid.dx);

% Plot the contour.
[ garbage, hCurve ] = contour(grid.xs{:}, theta - offset, [ 0 0 ], 'b-');

%---------------------------------------------------------------------------
% Add title if necessary.
if(nargin > 4)
  title(titleStr);
end

drawnow;
