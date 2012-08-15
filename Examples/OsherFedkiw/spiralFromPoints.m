function data = spiralFromPoints(grid, samples)
% spiralFromPoints: Generate an implicit surface function for a wound spiral.
%
%   data = spiralFromPoints(grid, samples)
%
% Generates an implicit surface function for a wound spiral.  This routine
%   is an attempt to duplicate the wound spiral interface found in Osher &
%   Sethian 88 and many subsequent level set publications discussing
%   curvature dependent flow; for example, O&F figure 4.1.
%
% The spiral is generated as a signed distance function to a set of
%   points from a parametric description of the spiral given in Osher &
%   Sethian 88:
%
%         \gamma(s) = (0.1 * exp(-10 * y(s)) - 0.05 * (0.1 - x(s)))
%                         * [ cos(a(s)); sin(a(s)) ],
%
%   where
%         a(s) = 25.0  * atan(10 * y(s))
%         x(s) =  0.1  * cos(2 * pi * s) + 0.1
%         y(s) =  0.05 * sin(2 * pi * s) + 0.1
%           s \in [ 0, 1 ]
%
% The sampling of the range of s mentioned in Osher & Sethian 88 is 
%   \delta s = 0.005 (200 points).
%
%
% Parameters:
%   grid         Grid on which the spiral will be defined.
%   samples      How many samples should be taken?
%
%   data         Implicit surface function for the wound spiral.
  
% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell 2/19/04

% This function cannot run until we have code to generate an implicit
%   surface function from a point cloud.
error('This function is not currently operational');
  
  if(grid.dim ~= 2)
    error('Spiral can only be created on two dimensional grids');
  end
  
  s = linspace(0, 1, samples)';
  x = 0.1  * cos(2 * pi * s) + 0.1;
  y = 0.05 * sin(2 * pi * s) + 0.1;
  a = 25 * atan(10 * y);
  
  multiplier = (0.1 * exp(-10 * y) - 0.05 * (0.1 - x));
  gamma =  repmat(multiplier, 1, 2) .* [ cos(a), sin(a) ];

  circle = 0.5 * [ cos(2 * pi * s), sin(2 * pi * s) ];
  inside = [ 0, 0 ];

  % Get a distance function.
  data = unsignedDistanceFromPoints(grid, circle);
  
  % Make it signed distance.
  data = flipSignInsideSurface(grid, data, inside);
