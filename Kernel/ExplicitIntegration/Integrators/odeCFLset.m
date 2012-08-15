function options = odeCFLset(varargin)
% odeCFLset Create/alter options for CFL constrained ode integrators.
%
%  options = odeCFLset('name1', value1, 'name2', value2, ...)
%  options = odeCFLset(oldopts, 'name1', value1, ...)
%
% Creates a new options structure (or alters an old one) for CFL
%   constrained ODE integrators.  Basically the same as Matlab's odeset
%   but with not nearly as many options.
%
% Without any options, it displays the available options, values and defaults.
%
% Available options (options names are case insensitive):
%
%   factorCFL    Scalar by which to multiply CFL timestep bound in order
%                  to determine the timestep to actually take.
%                  Typically in range (0,1), default = 0.5
%                  choose 0.9 for aggressive integration.
%
%   maxStep      Maximum step size (independent of CFL).
%                  Default is REALMAX.
%
%   postTimestep Function handle to a routine with prototype
%                       [ yOut schemeDataOut ] = f(t, yIn, schemeDataIn)
%                  which is called after every timestep and can be used
%                  to modify the state vector y or to modify or record
%                  information in the schemeData structure.
%                  Defaults to [], which calls no function.
%
%   singleStep   Specifies whether to exit integrator after a single
%                  CFL constrained timestep (for debugging).
%                  Either 'on' or 'off', default = 'off'.
%
%   stats        Specifies whether to display statistics.
%                  Either 'on' or 'off', default = 'off'.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/6/04

  %---------------------------------------------------------------------------
  % No output, no input means caller just wants a list of available options.
  if((nargin == 0) && (nargout == 0))
      fprintf('   factorCFL: [ positive scalar {0.5} ]\n');
      fprintf('     maxStep: [ positive scalar {REALMAX} ]\n');
      fprintf('postTimestep: [ function handle | {[]} ]\n');
      fprintf('  singleStep: [ on | {off} ]\n');
      fprintf('       stats: [ on | {off} ]\n');
      fprintf('\n');
      return;
  end
  
  %---------------------------------------------------------------------------
  % First input argument is an old options structure
  if((nargin > 0) && isstruct(varargin{1}))
    options = varargin{1};
    startArg = 2;
  else
    % Create the default options structure.
    options.factorCFL = 0.5;
    options.maxStep = realmax;
    options.postTimestep = [];
    options.singleStep = 'off';
    options.stats = 'off';
    startArg = 1;
  end
  
  %---------------------------------------------------------------------------
  % Loop through remaining name value pairs
  for i = startArg : 2 : nargin
    name = varargin{i};
    value = varargin{i+1};

    % Remember that the case labels are lower case.
    switch(lower(name))
     case 'factorcfl'
      if(isa(value, 'double') && (prod(size(value)) == 1) && ...
         (value > 0.0) && (value <= 1.0))
        options.factorCFL = value;
      else
        error('FactorCFL must be a scalar double value in range [ 0, 1 )');
      end
      
     case 'maxstep'
      if(isa(value, 'double') && (prod(size(value)) == 1) && (value > 0.0))
        options.maxStep = value;
      else
        error('MaxStep must be a positive scalar double value');
      end
            
     case 'posttimestep'
      if(isa(value,'function_handle'))
        options.postTimestep = value;
      else
        error('PostTimestep parameter must be a function handle');
      end

     case 'singlestep'
      if(isa(value, 'char') && (strcmp(value, 'on') || (strcmp(value, 'off'))))
        options.singleStep = value;
      else
        error('SingleStep must be one of the strings ''on'' or ''off''');
      end
      
     case 'stats'
      if(isa(value, 'char') && (strcmp(value, 'on') || (strcmp(value, 'off'))))
        options.stats = value;
      else
        error('Stats must be one of the strings ''on'' or ''off''');
      end
      
     otherwise
      error('Unknown odeCFL option %s', name);
      
    end
  end
