function data = signedDistanceIterative(grid, data0, tMax, errorMax, accuracy)
% signedDistanceIterative: Create a signed distance function iteratively.
%
%   data = signedDistanceIterative(grid, data0, tMax, errorMax, accuracy)
%
% Converts an implicit surface function into a signed distance function
%   by iterative solution of the reinitialization equation.
%
% Iterations continue to a fixed time or until the average relative change
%   in the function value between iterations drops low enough, whichever
%   comes first.
%
% In the reinitialization equation, information flows outward from the zero
%   level set at "speed" one, so to get a signed distance function in a band
%   of at least 10 grid cells around the zero level set, choose a minimum
%   tMax = 10 * max(grid.dx).
%
% Parameters:
%
%   grid	 Grid structure.
%   data0        Implicit surface function.
%   tMax         Time at which to halt the reinitialization iteration
%                  (default = norm(grid.max - grid.min)).
%   errorMax     If the average update of nodes drops below
%                  errorMax * max(grid.dx), then assume that 
%                  reinitialization has converged and return early 
%                  (default = 1e-3).
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%
%
%   data         Signed distance function.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/14/04

%---------------------------------------------------------------------------
% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% How aggressive should we be with the CFL condition?
factorCFL = 0.9;

%---------------------------------------------------------------------------
% Set defaults.
if(nargin < 3)
  tMax = max(grid.max - grid.min);
end

if(nargin < 4)
  errorMax = 1e-3;
end

if(nargin < 3)
  accuracy = 'medium';
end

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termReinit;
schemeData.grid = grid;
schemeData.initial = data0;

% Set up time approximation scheme.
%   Single step so that we can check convergence criterion.
integratorOptions = odeCFLset('factorCFL', factorCFL, 'singleStep', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  schemeData.derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  schemeData.derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  schemeData.derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  schemeData.derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end


%---------------------------------------------------------------------------
% Convergence criteria
deltaMax = errorMax * max(grid.dx) * prod(grid.N);

%---------------------------------------------------------------------------
% Reshape data array into column vector for ode solver call.
y = data0(:);

% Loop until tMax (subject to a little roundoff).
tNow = 0;
while(tMax - tNow > small * tMax)

  % Check for convergence (except for the first loop).
  if((tNow > 0) && (norm(y - y0, 1) < deltaMax))
    break;
  end
  
  % Always try to finish (in fact, single timesteps are taken).
  tSpan = [ tNow, tMax ];
  
  % Take a single timestep.
  y0 = y;
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

end

data = reshape(y, grid.shape);
