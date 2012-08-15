function [ ydot, stepBound ] = termRestrictSign(t, y, schemeData)
% termRestrictSign: restrict the sign of a term to be positive or negative.
%
% [ ydot, stepBound ] = termRestrictSign(t, y, schemeData)
%
% Given some other HJ term approximation, this function either restricts
%   the sign of that update to be either:
%       nonnegative (level set function only decreases), or
%       nonpositive (level set function only increases).
%
% The CFL restriction of the other HJ term approximation is returned 
%   without modification.
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
%   .innerFunc   Function handle for the approximation scheme to
%                  calculate the unrestricted HJ term.
%   .innerData   schemeData structure to pass to innerFunc.
%   .positive    Boolean, true if update must be positive (nonnegative)
%                  (optional, default = 1)
%
% It may contain addition fields at the user's discretion.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 3/26/04

  %---------------------------------------------------------------------------
  checkStructureFields(schemeData, 'innerFunc', 'innerData');

  %---------------------------------------------------------------------------
  % Default to positive (nonnegative) update restriction.
  if(isfield(schemeData, 'positive'))
    positive = schemeData.positive;
  else
    positive = 1;
  end

  %---------------------------------------------------------------------------
  % Get the unrestricted update.
  [ unRestricted, stepBound ] = ...
                      feval(schemeData.innerFunc, t, y, schemeData.innerData);

  %---------------------------------------------------------------------------
  % Restrict the update (stepBound is returned unchanged).  
  %   Do not negate for RHS of ODE (that is handled by innerFunc).
  if(positive)
    ydot = max(unRestricted, 0);
  else
    ydot = min(unRestricted, 0);
  end
