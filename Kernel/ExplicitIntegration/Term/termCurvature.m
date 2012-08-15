function [ ydot, stepBound ] = termCurvature(t, y, schemeData)
% termCurvature: approximate a motion by mean curvature term in an HJ PDE.
%
% [ ydot, stepBound ] = termCurvature(t, y, schemeData)
%
% Computes an approximation of motion by mean curvature for a
%   Hamilton-Jacobi PDE.  This is a second order equation that simplifies to
%   a heat equation if the function is a signed distance function.
%   Specifically:
%
%                 D_t \phi = - b(x) \kappa(x) \| \grad \phi \|.
%
%   where \kappa(x) is the mean curvature.
%
% Based on methods outlined in O&F, chapters 4.1 & 4.2.
%
% parameters:
%   t            Time at beginning of timestep.
%   y            Data array in vector form.
%   schemeData	 A structure (see below).
%
%   ydot	 Change in the data array, in vector form.
%   stepBound	 CFL bound on timestep for stability.
%
% schemeData is a structure containing data specific to this type of 
%   term approximation.  For this function it contains the field(s)
%
%   .grid	   Grid structure (see processGrid.m for details).
%   .curvatureFunc Function handle to finite difference curvature approx.
%                    Should provide both curvature and gradient magnitude.
%   .b		   Multiplier (should be non-negative for well-posedness);
%                    b may be a scalar constant or an array of size(data).
%
% It may contain addition fields at the user's discretion.
%
% schemeData.b can provide the multiplier in one of two ways:
%   1) For time invariant multipliers, a scalar or an array the same 
%      size as data.
%   2) For general multipliers, a function handle to a function with prototype
%      b = scalarGridFunc(t, data, schemeData), where the output b is the
%      scalar/array from (1) and the input arguments are the same as those
%      of this function (except that data = y has been reshaped to its
%      original size).  In this case, it may be useful to include additional
%      fields in schemeData.  
%
% In the notation of OF text,
%
%   data = y	  \phi
%   curvatureFunc function to calculate \kappa and |\grad \phi|
%   b		  b
%
%   delta = ydot  +b \kappa |\grad \phi|


% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell 6/3/03
% Calling parameters significantly modified, Ian Mitchell 2/13/04.

  %---------------------------------------------------------------------------
  checkStructureFields(schemeData, 'grid', 'curvatureFunc', 'b');

  %---------------------------------------------------------------------------
  grid = schemeData.grid;
  data = reshape(y, grid.shape);

  %---------------------------------------------------------------------------
  % Get multiplier
  if(isa(schemeData.b, 'double'))
    b = schemeData.b;
  elseif(isa(schemeData.b, 'function_handle'))
    b = feval(schemeData.b, t, data, schemeData);
  else
    error('schemeData.b must be a scalar, array or function handle');
  end

  %---------------------------------------------------------------------------
  % According to O&F equation (4.5).
  [ curvature, gradMag ] = feval(schemeData.curvatureFunc, grid, data);
  delta = -b .* curvature .* gradMag;
  
  %---------------------------------------------------------------------------
  % According to O&F equation (4.7).
  stepBound = 1 / (2 * max(b(:)) * sum(grid.dx .^ -2));
  
  % Reshape output into vector format and negate for RHS of ODE.
  ydot = -delta(:);
