function [ data, g, data0 ] = ...
                       nonconvexLF(accuracy, dissType, gridDim, gridSize, tMax)
% nonconvexLF: demonstrate Lax-Friedrichs on a nonconvex Hamiltonian.
%
%   [ data, g, data0 ] = nonconvexLF(accuracy, dissType, gridDim, gridSize, tMax)
%  
% This function demonstrates how the Lax-Friedrichs HJ term approximation
%   termLaxFriedrichs can be used to approximate the solution to a multiple
%   dimensional version of a general HJ PDE with an H function nonconvex in
%   \grad \phi.
%
% The example is taken from Osher & Shu, "High-Order Essentially
%   Nonoscillatory Schemes for Hamilton-Jacobi Equations",
%   SIAM J. Numerical Analysis, vol 28, num 4, pp 907-922 (August, 1991).
%   In particular, we are doing the nonconvex versions of examples 1 and 2.
%
% Specifically, in dimension d
%   PDE:  D_t \phi + H(x, D_x \phi) = 0
%   BC:   periodic in domain [ -d, +d ]
%   IC:   \phi(x,0) = -\cos((\pi / d) * \sum_i x_i)
%   Ham:  H(x,p) = -\cos(\alpha + \sum_i p_i)
%
% Parameters (all inputs have defaults):
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst (default).
%                  'medium'      Use odeCFL2 and upwindFirstENO2.
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   dissType     Which type of dissipation to use with the LF approximation?
%                  'global'      Use artificialDissipationGLF (default).
%                  'local'       Use artificialDissipationLLF.
%                  'locallocal'  Use artificialDissipationLLLF.
%   gridDim      Number of dimensions.  Default = 2.
%   gridSize     Number of nodes in each dimension. 
%                  Default = 80 / 2^(gridDim-1).
%   tMax         Final time.  Default = 1.5 / pi^2.
%
%   data         Implicit surface function at tMax.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at initial time.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/12/04

%---------------------------------------------------------------------------
% You will see some executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Set defaults.
if(nargin < 1)
  accuracy = 'low';
end

if(nargin < 2)
  dissType = 'global';
end

if(nargin < 3)
  gridDim = 2;
end

if(nargin < 4)
  gridSize = 80 / 2^(gridDim - 1);
end

if(nargin < 5)
  tMax = 1.5 / pi^2;
end

%---------------------------------------------------------------------------
% Integration parameters.
t0 = 0;                      % Start time.
alpha = 1.0;                 % Parameter in the Hamiltonian.

% Plotting choices.
plotIntermediate = 1;        % Plot at intermediate timesteps?
tPlot = tMax / 10;           % Period at which plot should be produced.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

%---------------------------------------------------------------------------
% Use periodic boundary conditions?
periodic = 1;

% Create the grid.
g.dim = gridDim;
g.min = -g.dim;
g.N = gridSize;
if(periodic)
  g.max = (gridSize - 2) / gridSize * g.dim;
  g.bdry = @addGhostPeriodic;
else
  g.max = g.dim;
  g.bdry = @addGhostExtrapolate;
end
g = processGrid(g);

%---------------------------------------------------------------------------
% What kind of display?
switch(g.dim)
 case 1
  displayType = 'plot';
 case 2
  displayType = 'surf';
 case 3
  displayType = 'slice';
 otherwise
  error('Default display type undefined for dimension %d', g.dim);
end

%---------------------------------------------------------------------------
% Create initial conditions.
%   Extending the 1 and 2 dimensional examples to higher dimension:
%       \phi(x, 0) = -\cos((\pi / dimension) * \sum_i x_i)
data = g.xs{1};
for i = 2 : g.dim
  data = data + g.xs{i};
end
data = -cos((pi / g.dim) * data);
data0 = data;

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.grid = g;
schemeData.hamFunc = @nonconvexHamFunc;
schemeData.partialFunc = @nonconvexPartialFunc;

% This parameter is not required by termLaxFriedrichs, 
%   but is used by hamFunc and partialFunc.
schemeData.alpha = alpha;

%---------------------------------------------------------------------------
% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
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
% What kind of dissipation?
switch(dissType)
 case 'global'
  schemeData.dissFunc = @artificialDissipationGLF;
 case 'local'
  schemeData.dissFunc = @artificialDissipationLLF;
 case 'locallocal'
  schemeData.dissFunc = @artificialDissipationLLLF;
 otherwise
  error('Unknown dissipation function %s', dissFunc);
end
  
%---------------------------------------------------------------------------
% Initialize Display
f = figure;

h = visualizeLevelSet(g, data, displayType, level, ['t = ' num2str(t0)]);

hold on;
if(g.dim > 1)
  axis(g.axis);
  daspect([ 1 1 1 ]);
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  if(plotIntermediate)
    if(pauseAfterPlot)
      % Wait for last plot to be digested.
      pause;
    end
    
    % Get correct figure.
    figure(f);

    % Delete last visualization if necessary.
    if(deleteLastPlot)
      delete(h);
    end

    % Create new visualization.
    h = visualizeLevelSet(g, data, displayType, level, ['t = ' num2str(tNow)]);

  end
end

endTime = cputime;
fprintf('Total execution time %g seconds\n', endTime - startTime);

%---------------------------------------------------------------------------
% Plot the final result, if it has not already appeared.
if(~plotIntermediate)
  % Get correct figure.
  figure(f);
  
  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end
  
  % Create new visualization.
  visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
end
  


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = nonconvexHamFunc(~, data, deriv, schemeData)
% nonconvexHamFunc: analytic Hamiltonian for nonconvex HJ PDE.
%
% hamValue = nonconvexHamFunc(t, data, deriv, schemeData)
%
% This function implements the hamFunc prototype to calculate the analytic
%   Hamiltonian for the nonconvex case from Examples 1 & 2 of Osher & Shu
%   1991.
%
% Specifically, H(x,p) = -\cos(\alpha + \sum_i p_i)
%
% Parameters:
%   t            Time at beginning of timestep (ignored).
%   data         Data array.
%   deriv	 Cell vector of the costate (\grad \phi).
%   schemeData	 A structure (see below).
%
%   hamValue	 The analytic hamiltonian.
%
% schemeData is a structure containing data specific to this Hamiltonian
%   For this function it contains the field(s):
%
%   .grid	 Grid structure.
%   .alpha       Parameter in the Hamiltonian.
%
% Ian Mitchell 2/12/04

checkStructureFields(schemeData, 'grid', 'alpha');

hamValue = schemeData.alpha * ones(size(data));
for i = 1 : schemeData.grid.dim;
  hamValue = hamValue + deriv{i};
end
hamValue = -cos(hamValue);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = ...
             nonconvexPartialFunc(~, ~, ~, ~, schemeData, ~)
% nonconvexPartialFunc: Hamiltonian partial function for nonconvex example.
%
% alpha = nonconvexPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
%
% This function implements the partialFunc prototype to calculate the
%   extrema of the absolute value of the partials of the analytic
%   Hamiltonian with respect to the costate (gradient) for the nonconvex
%   case from Examples 1 & 2 of Osher & Shu 1991.
%
% NOTE: There are two "alpha" in this function: the return value alpha
%   and the Hamiltonian parameter schemeData.alpha.  Do not mix them up!
%
% Specifically, max_p | \partial H(x,p) / \partial p_dim |,
%   where H(x,p) = -\cos(\alpha + sum_i p_i).
%
% So we are maximizing over p the function |\sin(\alpha + \sum_i p_i)|.
%   Ideally, we would solve this nonconvex optimization.
%   For simplicity, we'll just use |sin(...)| <= 1.
%
% Parameters:
%   t            Time at beginning of timestep (ignored).
%   data         Data array.
%   derivMin	 Cell vector of minimum values of the costate (\grad \phi).
%   derivMax	 Cell vector of maximum values of the costate (\grad \phi).
%   schemeData	 A structure (see below).
%   dim          Dimension in which the partial derivatives is taken.
%
%   alpha	 Maximum absolute value of the partial of the Hamiltonian
%		   with respect to the costate in dimension dim for the 
%                  specified range of costate values (O&F equation 5.12).
%		   Note that alpha can (and should) be evaluated separately
%		   at each node of the grid.
%
% schemeData is a structure containing data specific to this Hamiltonian
%   For this function it contains the field(s):
%
%   .alpha       Parameter in the Hamiltonian.
%
% Ian Mitchell 2/12/04

checkStructureFields(schemeData, 'alpha');

% We want to maximize |\sin(...)| over a range of p;
%   conservative choice is just one.

alpha = 1;
