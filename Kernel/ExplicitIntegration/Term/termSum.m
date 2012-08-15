function [ ydot, stepBound ] = termSum(t, y, schemeData)
% termSum: Combine a collection of spatial HJ term approximations.
%
% [ ydot, stepBound ] = termSum(t, y, schemeData)
%
% This function independently evaluates a collection of HJ term 
%   approximations and returns their elementwise sum.
%
% The CFL restrictions are inverse summed:
%	stepBound_sum = (\sum_i (1 / stepBound_i))^-1
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
%   .innerFunc   A cell vector of function handles to the term 
%                  approximations that will be summed.
%   .innerData   A cell vector of schemeData structures to pass to 
%                  each of the innerFunc elements.  This vector
%                  must be the same size as the vector innerFunc.
%
% It may contain addition fields at the user's discretion.
%

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 4/9/04

  %---------------------------------------------------------------------------
  checkStructureFields(schemeData, 'innerFunc', 'innerData');

  %---------------------------------------------------------------------------
  % Check that innerFunc and innerData are the same size cell vectors.
  if(~iscell(schemeData.innerFunc) | ~iscell(schemeData.innerData))
    error('schemeData.innerFunc and schemeData.innerData %s', ...
          'must be cell vectors');
  end

  numSchemes = length(schemeData.innerFunc(:));

  if(numSchemes ~= length(schemeData.innerData(:)))
    error('schemeData.innerFunc and schemeData.innerData must be %s', ...
          'the same length');
  end

  %---------------------------------------------------------------------------
  % Calculate sum of updates (inverse sum of stepBounds).
  ydot = 0;
  stepBoundInv = 0;
  for i = 1 : numSchemes
    [ updateI stepBoundI ] = feval(schemeData.innerFunc{i}, t, y, ...
                                   schemeData.innerData{i});
    ydot = ydot + updateI;
    stepBoundInv = stepBoundInv + 1 / stepBoundI;
  end
  stepBound = 1 / stepBoundInv;
