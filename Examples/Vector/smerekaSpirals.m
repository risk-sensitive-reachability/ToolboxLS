function [ dataCurve, dataMask, g, tPlot ] = ...
                            smerekaSpirals(whichFig, exactCopy, accuracy, tMax)
% smerekaSpirals: example of dynamic open curves by vector level sets.
%
%   [ dataCurve, dataMask, g, tPlot ] = ...
%                           smerekaSpirals(whichFig, exactCopy, accuracy, tMax)
%  
% This function demonstrates use of vector level set methods for dynamic
% implicit open curves by recreating examples of spiral crystal growth from
%
%      Peter Smereka, "Spiral Crystal Growth," Physica D 138, pp. 282-301
%      (2000).
%
% Note that this function may take a LONG TIME to run because the motion by
% mean curvature dictates a very small timestep.  In 2007, it took more
% than 24 hours to complete the co-rotating spirals.
%
% The essential idea for open curves is to propagate two implicit surface
% functions at once.  One tracks a closed curve and the other a masking
% region.  Where the mask function is positive, the curve is considered
% real, and where the mask function is negative the curve is a dummy. The
% ends of the open real curve lie where the mask function is zero.
%
% The dynamics for curve motion are essentially motion in the normal
% direction plus motion by curvature.
%
% This function was originally designed as a script file, so most of the
% options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
% flow field parameters, etc.
%
% Input Parameters (all inputs have defaults):
%
%   whichFig: Number specifying which figure from Smereka's paper to
%   recreate.  Default is 4.
%                    4           Counter-rotating spirals, large spacing.
%                    7           Counter-rotating spirals, small spacing.
%                    12          Co-rotating spirals.
%                    17          Single spiral.
%
%   exactCopy: Boolean indicating whether Smereka's exact figure should be
%   recreated.  If not, use the same parameters and initial conditions, but
%   plot equally spaced time increments.  Default is 0.
%
%   accuracy: Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%
%   tMax: Final time.  Data arrays are returned at this time. In order to
%   see a full figure recreation, this value must be greater than the time
%   of the final corresponding subplot from Smereka's paper.  Default is
%   the time of the corresponding final subplot.
%
% Output parameters:
%
%   dataCurve: Implicit surface function for the curve. Corresponds to
%   Smereka's \phi(x,t). Optionally (depending on whether output tPlot is
%   requested) returns a cell vector, each element of which is the implicit
%   surface function at a sample time corresponding with one element of
%   tPlot.
%
%   dataMask: Implicit surface function for the masking set. Corresponds to
%   Smereka's \psi(x,t). Optionally (depending on whether output tPlot is
%   requested) returns a cell vector, each element of which is the implicit
%   surface function at a sample time corresponding with one element of
%   tPlot.
%
%   g: Grid structure on which data was computed.
%
%   tPlot: Column vector of sample times at which plots were produced.
%   Requesting this output indicates that dataCurve and dataMask should
%   return their cell vector versions.

% Copyright 2005 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/08/05

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
g.N = 401;
g.min = -10 * ones(g.dim, 1);
g.max = +10 * ones(g.dim, 1);
g.bdry = @addGhostExtrapolate;
%g.bdry = @addGhostNeumann;
g = processGrid(g);

%---------------------------------------------------------------------------
% Problem parameters.

% Mollification parameter (for sign function), suggested by Smereka.
epsilon = 4 * max(g.dx);

% By default, recreate figure 10.5 reprinted in Osher & Fedkiw text.
if(nargin < 1)
  whichFig = 4;
end

% Figure parameters from Smereka.
switch(whichFig)

 case 4,
  lambda = 0.1;
  tLast = 11;
  tPlot = [ 1, 5, 7, 11 ];
  curve0 = g.xs{2};
  mask0 = 4 - abs(g.xs{1});

 case 7,
  lambda = 0.2;
  tLast = 25;
  tPlot = [ 1, 8, 16, 25 ];
  curve0 = g.xs{2};
  mask0 = 1 - abs(g.xs{1});

 case 12,
  lambda = 0.2;
  tLast = 25;
  tPlot = [ 1, 7, 16, 25 ];
  curve0 = -min(1 + g.xs{2}, 1 - g.xs{2});
  % The equation (21) for this function was buggy in Smereka.
  mask0 = max((g.xs{1} - 10)/10 + g.xs{2}, (10 - g.xs{1})/10 - g.xs{2} - 2);

 case 17,
  lambda = 0.2;
  tLast = 25;
  tPlot = [ 1, 7, 16, 25 ];
  curve0 = -g.xs{2};
  mask0 = g.xs{1};

 otherwise,
  error([ 'Unknown whichFig value: ' num2str(whichFig) ]);

end

data0 = { curve0; mask0 };
data = data0;

%---------------------------------------------------------------------------
% Integration parameters.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% How many steps of reinitialization after each timestep?
%   Set to zero for no reinitialization (which seems to be insufficient).
reinitSteps = 1;

% Should we show just the open curves,
%   or both the mask set and the entire contour?
showAll = 0;

%---------------------------------------------------------------------------
% We are not viewing level sets, so choose level to disable warnings.
level = [];

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 0;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

%---------------------------------------------------------------------------
% Decide whether we are recreating the exact Smereka plots, or
%   just using their initial conditions and parameters.
if(nargin < 2)
  exactCopy = 0;
end

if(nargin < 4)
  tMax = tLast;
end

if(exactCopy)
  % Adjust tPlot vector to take account of tMax.
  while(tMax < tPlot(end))
    tPlot = tPlot(1:end-1);
  end
  if(tMax > tPlot(end))
    tPlot = [ tPlot, tMax ];
  end

  plotSteps = length(tPlot);
  deleteLastPlot = 0;
  useSubplots = 1;

else
  tPlot = linspace(t0, tMax, plotSteps);

end
  
% What version of the output is requested?
if(nargout > 3)
  keepPlotData = 1;
  dataCurve = cell(plotSteps, 1);
  dataMask = cell(plotSteps, 1);
else
  keepPlotData = 0;
end

% Indexing constants for the two level set functions.
iCurve = 1;
iMask = 2;

% How many level set functions?
iCount = 2;

%---------------------------------------------------------------------------
if(nargin < 3)
  accuracy = 'medium';
end

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

%---------------------------------------------------------------------------
% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.8, 'stats', 'on');

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Motion is exactly symmetric:
%   Both level set functions move under the same dynamics and
%   interact symmetrically with one another.
%   Consequently, we need only one schemeFunc and one schemeData.
%
% For more general motions, either or both can be cell vectors.

% Motion for both sets will consist of a sum of other terms.
schemeFunc = @termSum;

% Motion in the normal direction at speed 1 with sign set
%   by the other level set function.
% Since the latter changes with time, we need a time dependent speed function.
normalData.grid = g;
normalData.derivFunc = derivFunc;
normalData.speed = @maskedSpeedFunc;
normalData.maskFunc = @smoothSign;
normalData.speedConst = +1;
normalData.epsilon = epsilon;
normalData.passVLS = 1;

% Motion by curvature at speed lambda.
%   Note that Smereka's paper uses a rather confusing formulation 
%   (eg (9) & (10)) that results in the curvature speed being
%   \lambda * \sign(\phi)^2 \neq \lambda (for the smoothed sign function).
curveData.grid = g;
curveData.curvatureFunc = @curvatureSecond;
curveData.b = @maskedSpeedFunc;
curveData.maskFunc = @smoothSignSquared;
curveData.speedConst = lambda;
curveData.epsilon = epsilon;
curveData.passVLS = 1;

% All motion consists of the sum of the same two terms.
schemeData = cell(iCount, 1);
schemeData{iCurve}.innerFunc = { @termNormal; @termCurvature };
schemeData{iCurve}.innerData = { normalData; curveData };

% Only difference for the masking function is that the normal motion is
%   in the opposite direction (despite what Smereka (10) might say).
normalData.speedConst = -normalData.speedConst;
schemeData{iMask}.innerFunc = { @termNormal; @termCurvature };
schemeData{iMask}.innerData = { normalData; curveData };

% We also need access to the grid structure at the termSum level
%   so that termNormal can reshape the arrays properly.
for i = 1 : iCount
  schemeData{i}.grid = g;
end

%---------------------------------------------------------------------------
% Set up reinitialization if necessary.
if(reinitSteps)
  % Register reinitialization postTimestep routine.
  integratorOptions = odeCFLset(integratorOptions, ...
                                'PostTimestep', @postTimestepReinit);
  
  % Add necessary components to schemeData structures.
  for i = 1 : iCount
    schemeData{i}.reinitAccuracy = accuracy;
    schemeData{i}.reinitSteps = reinitSteps;
  end
  
end

%---------------------------------------------------------------------------
% Initialize Display
f = figure;
plotCount = 1;

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  subplot(rows, cols, plotCount);
end

% Is the initial time a plotting point?
if(tPlot(1) == t0)
  [ hCurve, hMask ] = visualizeOpenCurve(g, data{iCurve}, data{iMask}, ...
                                         showAll, [ 't = ' num2str(t0) ]);
  if(keepPlotData)
    dataCurve{plotCount} = data{iCurve};
    dataMask{plotCount} = data{iMask};
  end
  plotCount = plotCount + 1;
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
y0 = cell(iCount, 1);
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  for i = 1 : iCount
    y0{i} = data{i}(:);
  end

  % How far to step?
  tSpan = [ tNow, tPlot(plotCount) ];
  
  % Take a timestep.
  [ t, y, schemeData ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                               integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  for i = 1 : iCount
    data{i} = reshape(y{i}, g.shape);
  end

  % Record plot data if necessary.
  if(keepPlotData)
    dataCurve{plotCount} = data{iCurve};
    dataMask{plotCount} = data{iMask};
  end

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % Create new visualizations.
  figure(f);
  if(useSubplots)
    subplot(rows, cols, plotCount);
  end

  if(deleteLastPlot)
    delete(hCurve);
    delete(hMask);
  end
  [ hCurve, hMask ] = visualizeOpenCurve(g, data{iCurve}, data{iMask}, ...
                                         showAll, [ 't = ' num2str(tNow) ]); 
  hold on;
  plotCount = plotCount + 1;

end

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);

if(~keepPlotData)
  % Data for individual plots was not maintained; return only the final data.
  dataCurve = data{iCurve};
  dataMask = data{iMask};
end



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function speed = maskedSpeedFunc(t, data, schemeData)
% maskedSpeedFunc: Constant speed, direction dictated by auxiliary function.
%
%   speed = maskedSpeedFunc(t, data, schemeData)
%
% Implements the scalarGridFunc prototype to provide a speed function for
%   some kind of motion.  In this case, the speed is masked against some
%   function which depends on the second element of the vector level set
%   being evolved.
%
% This type of motion appears in the evolution of open surfaces in 
%   Smereka, (9) and (10).
%
% It is assumed that we are determining the motion of the implicit surface
%   function stored in data{1} using data from the implicit surface
%   function stored in data{2}.  Note that data{1} and data{2} must be 
%   the same size array.
%
% Parameters:
%   t            Current time.
%   data         Vector level set function (a cell vector).
%   schemeData   Structure or cell vector of structures (see below).
%
%   speed        Array containing the speed of motion at each grid node.
%
% schemeData (or schemeData{1} if schemeData is a cell vector)
%   is a structure containing data specific to this type of 
%   term approximation.  For this function it contains the field(s)
%
%   .speedConst  Basic speed of the motion.  May be scalar or an
%                  array the same size as data{1}.
%   .maskFunc    Function handle for the function which generates
%                  the result of masking.
%   .epsilon     Smoothing factor for the mask function.
%
% schemeData may contain other fields.

  % Note that we are receiving a cell vector of data and (possibly) schemeData.
  %   The first element of these cell vectors is the active one.
  % We don't actually need the other elements of schemeData, but we
  %   do need the other elements of data, so we have no choice.

  if(iscell(schemeData))
    thisSchemeData = schemeData{1};
  else
    thisSchemeData = schemeData;
  end

  checkStructureFields(thisSchemeData, 'speedConst', 'maskFunc', 'epsilon');

  speed = thisSchemeData.speedConst .* ...
          feval(thisSchemeData.maskFunc, data{2}, thisSchemeData.epsilon);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function sgn = smoothSign(phi, epsilon)
% smoothSign: a smoothed version of the sign function.
%
%   sgn = smoothSign(phi, epsilon)
%
% This function computes the smoothed version of the sign function
%   suggested in Smereka (16).
%
%                                   +1,           \phi  >    \epsilon
%        \sgn_\epsilon(\phi) = \phi / \epsilon,  |\phi| \leq \epsilon 
%                                   -1,           \phi  <   -\epsilon
% Parameters:
%   phi          The argument of the sign function.
%   epsilon      The smoothing factor.
%
%   sgn          The smoothed sign of the argument \sgn_\epsilon(\phi).

sgn = min(max(phi / epsilon, -1), +1);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function sgn2 = smoothSignSquared(phi, epsilon)
% smoothSignSquared: a smoothed version of the sign function, squared.
%
%   sgn2 = smoothSignSquared(phi, epsilon)
%
% This function computes the square of the smoothed version of the sign 
%   function suggested by Smereka.  Note that this function is not
%   the identity because of the smoothing.
%
% Parameters:
%   phi          The argument of the sign function.
%   epsilon      The smoothing factor.
%
%   sgn2         The square of the smoothed sign of the argument.

sgn2 = smoothSign(phi, epsilon).^2;
