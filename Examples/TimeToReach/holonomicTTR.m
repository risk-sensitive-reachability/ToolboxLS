function [ mttr, attr, data, gridOut, data0 ] = ...
                                       holonomicTTR(whichNorm, accuracy,gridIn)
% holonomicTTR: demonstrate a holonomic time to reach function.
%
%   [ mttr, attr, data, gridOut, data0 ] = ...
%                                     holonomicTTR(whichNorm, accuracy, gridIn)
%  
% In this example we calculate the minimum time to reach a small ball
%   around the origin when motion is allowed in any direction.
%
% The dynamics are
%
%        \dot x    = b_1
%	 \dot y    = b_2
%
%   where input ||b|| \leq 1 is trying to hit the target.  The norm
%   in which b is bounded can be modified to produce different shapes
%   of time to reach function.
%
% The goal is slightly more ambitious than the standard reach set, since
%   we would also like to record the minimum time to reach function by
%   tracking the time at which the reach set arrives at each node.
%
% Reinitialization is available, but only seems to degrade solution quality.
%
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   aircraft parameters, etc.
%
% Parameters:
%
%   whichNorm    Controls the norm in which b is bounded.  Default is 2.
%                  '1',1,'sum'       norm(b, 1) = sum(abs(b)) <= 1.
%                  '2',2,'rms'       norm(b, 2) = sqrt(sum(b.^2)) <= 1.
%                  'inf',inf,'max'   norm(b, inf) = max(abs(b)) <= 1.
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   gridIn       Computational grid.
%                  Default is dim = 2; min = -1; max = +1; N = 101;
%                  and bdry = @addGhostExtrapolate.
%
%   mttr         Minimum time to reach function at t_max.
%   attr         Analytic minimum time to reach function.
%                  Determined by calling analyticHolonomic 
%                  or analyticSumSquare.
%   data         Implicit surface function at t_max for the reach set.
%   gridOut      Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 12/06/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.

%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.2;                  % End time.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% What kind of dissipation?
dissType = 'global';

%---------------------------------------------------------------------------
% Problem Parameters.

if(nargin < 1)
  whichNorm = 2;
end

% Bounds on the input magnitude.
inputBound = 1;

% Radius of the target set.
%   A value of zero will be translated into a value of (smallTarget * g.dx).
targetRadius = 0.1;

% What is considered a small target?
smallTarget = 1.0;

% Should we reinitialize after every step?
reinitialize = 0;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Visualize the reachable set.
displayType = 'contour';

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;

% Show analytic solution for comparison?
showAnalytic = 1;

% What type of initial conditions?
whichIC = 'circle';
whichIC = 'square';
%whichIC = 'analytic';

%---------------------------------------------------------------------------
% Create the grid (if necessary).

if(nargin < 3)
  g.dim = 2;
  g.min = -1;
  g.max = +1;
  g.bdry = @addGhostExtrapolate;
  g.N = 101;
else
  g = gridIn;
end

g = processGrid(g);

if(nargout > 3)
  gridOut = g;
end

% Check the target radius
if(targetRadius == 0)
  % Precisely zero target radius will not work.  
  %   Make it small relative to the size of the grid.
  targetRadius = smallTarget * max(g.dx);
end

%---------------------------------------------------------------------------
% Evaluate analytic solution for target of appropriate radius.
attr = analyticHolonomicTTR(whichNorm, g) - targetRadius;

% Alternative.  Should only be used with 'square' IC option and 'sum' norm.
attr = analyticSumSquareTTR(targetRadius, g);

% Create initial conditions.
switch(whichIC)
 case 'circle'
  data = shapeSphere(g, zeros(g.dim, 1), targetRadius);

 case 'square'
  data = shapeRectangleByCenter(g, zeros(g.dim, 1), 2 * targetRadius);

 case 'analytic'
  data = attr;

 otherwise
  error('Unknown type of initial conditions whichIC = ', whichIC);
end

if(nargout > 4)
  data0 = data;
end

% Although we needed negative values in attr to properly define the
%   implicit surface for the 'analytic' IC case, real time to reach
%   functions do not have negative values.
attr = max(0, attr);

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.hamFunc = @holonomicHamFunc;
schemeData.partialFunc = @holonomicPartialFunc;
schemeData.grid = g;

% The Hamiltonian and the partials need problem parameters.
schemeData.inputBound = inputBound;
schemeData.whichNorm = whichNorm;

%---------------------------------------------------------------------------
% Choose degree of dissipation.

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
if(nargin < 2)
  accuracy = 'medium';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.75, 'stats', 'on');

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

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Restrict the Hamiltonian so that reachable set only grows.
%   The Lax-Friedrichs approximation scheme MUST already be completely set up.
innerFunc = schemeFunc;
innerData = schemeData;
clear schemeFunc schemeData;

% Wrap the true Hamiltonian inside the term approximation restriction routine.
schemeFunc = @termRestrictUpdate;
schemeData.innerFunc = innerFunc;
schemeData.innerData = innerData;
schemeData.positive = 0;

%---------------------------------------------------------------------------
% Set up minimum time to reach recording using postTimestepFunc.
integratorOptions = odeCFLset(integratorOptions, ...
                              'postTimestep', @postTimestepTTR);

% Initialize the minimum time to reach function by calling
%   the postTimestepFunc once with initial data.
y = data(:);
[ y, schemeData ] = feval(@postTimestepTTR, t0, y, schemeData);
data = reshape(y, g.shape);

%---------------------------------------------------------------------------
if(reinitialize)

  % We need some more parameters to drive the reinitialization.
  %   But most can just be the defaults.
  schemeData.grid = g;
  schemeData.reinitAccuracy = accuracy;
  schemeData.reinitSteps = 1;

  % Register the entire cell vector of postTimestep routines.
  integratorOptions = odeCFLset(integratorOptions, 'postTimestep', ...
                                { @postTimestepReinit; @postTimestepTTR });
end

%---------------------------------------------------------------------------
% Initialize Display
f = figure;

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);
hold on;

if(showAnalytic)
  ha = [];
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
  [ t, y, schemeData ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                               integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % Get correct figure, and remember its current view.
  figure(f);

  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end

  % Move to next subplot if necessary.
  if(useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  if(showAnalytic)
    delete(ha);
    ha = visualizeLevelSet(g, attr, displayType, tNow);
    set(ha, 'LineStyle', '--');
  end

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);

end

endTime = cputime;
fprintf('Total execution time %g seconds\n', endTime - startTime);

% Extract the minimum time to reach function from the schemeData structure.
%   Reshape it into an array, and replace the +inf entries with NaN to
%   make the visualization more pleasant.
mttr = reshape(schemeData.ttr, g.shape);
unreached = find(isinf(mttr));
mttr(unreached) = NaN;



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = holonomicHamFunc(t, data, deriv, schemeData)
% holonomicHamFunc: analytic Hamiltonian for the holonomic test case.
%
% hamValue = holonomicHamFunc(t, data, deriv, schemeData)
%
% This function implements the hamFunc prototype for the holonomic test case:
%
%   hamValue = -(input1 * deriv{1} + input2 * deriv{2})
%
% The actual values of input1 and input2 depend on deriv and on the
%   type of norm in which the input is bounded.
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
%   .inputBound  Bound on the magnitude of the input.
%   .whichNorm   Norm in which the input is bounded.
%
% Ian Mitchell 12/06/04

checkStructureFields(schemeData, 'inputBound', 'whichNorm');

switch(schemeData.whichNorm)

 case { '1'; 1; 'sum' }
  hamValue = -(-schemeData.inputBound * max(abs(deriv{1}), abs(deriv{2})));

 case { '2'; 2; 'rms' }
  hamValue = -(-schemeData.inputBound * sqrt(deriv{1} .^2 + deriv{2} .^2));

 case { 'inf'; inf; 'max' }
  hamValue = -(-schemeData.inputBound * abs(deriv{1}) - ...
                schemeData.inputBound * abs(deriv{2}));

 otherwise
  error('Unknown type of norm: %s', schemeData.whichNorm);
end



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = ...
      holonomicPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
% holonomicPartialFunc: partial fcn for the holonomic test case.
%
%   alpha = holonomicPartialFunc(t, data, derivMin, ...
%                                       derivMax, schemeData, dim)
%
% This function implements the partialFunc prototype for 
%   the holonomic test case.
%
% It calculates the extrema of the absolute value of the partials of the 
%   analytic Hamiltonian with respect to the costate (gradient).
%
% Referring to the Hamiltonian from the previous function, 
%   a loose bound on those partials is:
%
%   in dimension 1: schemeData.inputBound
%   in dimension 2: schemeData.inputBound
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
%   .inputBound  Bound on the magnitude of the input.
%
% Ian Mitchell 12/06/04

checkStructureFields(schemeData, 'inputBound');
alpha = schemeData.inputBound;
