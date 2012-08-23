function [ data, g, data0 ] = maskDemo(accuracy, displayType)
% maskDemo: demonstrate the masking process on a convective flow.
%
% [ data, g, data0 ] = maskDemo(accuracy, displayType)
%  
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters:
%
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst (default).
%                  'medium'      Use odeCFL2 and upwindFirstENO2.
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   displayType  String to specify how to display results.
%                  The specific string depends on the grid dimension;
%                  look at the helper visualizeLevelSet to see the options
%                  (optional, default depends on grid dimension).
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 4/11/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Which of the tasks do you wish to perform?
doMask = 1;
doMin = 1;

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.0;                  % End time.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 0;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

%---------------------------------------------------------------------------
% Use periodic boundary conditions?
periodic = 0;

% Create the grid.
g.dim = 2;
g.min = -1;
g.dx = 1 / 50;
if(periodic)
  g.max = (1 - g.dx);
  g.bdry = @addGhostPeriodic;
else
  g.max = +1;
  g.bdry = @addGhostExtrapolate;
end
g = processGrid(g);

%---------------------------------------------------------------------------
% What kind of display?
if(nargin < 2)
  switch(g.dim)
   case 1
    displayType = 'plot';
   case 2
    displayType = 'contour';    
   case 3
    displayType = 'surface';
   otherwise
    error('Default display type undefined for dimension %d', g.dim);
  end
end

%---------------------------------------------------------------------------
% Create initial conditions (a circle/sphere).
%   Note that in the periodic BC case, these initial conditions will not
%   be continuous across the boundary unless the circle is perfectly centered.
%   In practice, we'll just ignore that little detail.
center = [ -0.4; 0.0; 0.0; 0.0 ];
radius = 0.35;
data = shapeSphere(g, center, radius);

%---------------------------------------------------------------------------
% Create mask (a smaller circle/sphere).
maskCenter = [ 0.0; 0.0; 0.0; 0.0 ];
maskRadius = 0.15;
mask = shapeSphere(g, maskCenter, maskRadius);

% The moving set can be anywhere outside the masked region.
%   So what we really need is the complement of the masked region.
mask = -mask;

% Need to ensure that the initial conditions satisfy the mask.
data = max(data, mask);

data0 = data;

%---------------------------------------------------------------------------
% Most of the time in constant flow case, we want flow in a 
%   distinguished direction, so assign first dimension's flow separately.
v = 0 * ones(g.dim);
v(1) = 1.5;
v = num2cell(v);

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'medium';
end

% Set up spatial term approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = v;
schemeData.grid = g;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

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
% Set up data required for the mask operation.
%   Mask will be compared to vector form of data array used by integrator.
schemeData.mask = mask(:);
schemeData.doMask = doMask;

% Also keep track of minimum of phi over time.
%   Minimum will be kept in vector form used by integrator.
schemeData.min = data(:);
schemeData.doMin = doMin;

% Let the integrator know what function to call.
integratorOptions = odeCFLset(integratorOptions, ...
                              'postTimestep', @maskAndKeepMin);

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
  %   Record returned schemeData structure to keep track of min over time.
  [ t y schemeData ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
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

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);

end

%---------------------------------------------------------------------------
endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);

%---------------------------------------------------------------------------
% Process and display final results.
minOverTime = reshape(schemeData.min, g.shape);

if(g.dim == 2)
  % Display initial set, mask, minimum over time, and final set.
  figure;
  lev = [ level level ];
  [ ~, hI ] = contour(g.xs{1}, g.xs{2}, data0, lev, 'b--');
  hold on;
  [ ~, hF ] = contour(g.xs{1}, g.xs{2}, data, lev, 'r-');
  [ ~, hT ] = contour(g.xs{1}, g.xs{2}, minOverTime, lev, 'k:');
  [ ~, hM ] = contour(g.xs{1}, g.xs{2}, mask, lev, 'g-.');

  hs = [ hI; hF; hT; hM ];

  set(hs, 'LineWidth', 2.0);

  legend([ hI(1), hF(1), hT(1), hM(1) ], ...
         {'initial', 'final', 'min over time', 'mask'}, 2);
  axis equal
  axis(g.axis);
else
  warning('Cannot create final plot in dimensions other than 2D'); %#ok<WNTAG>
end



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function [ yOut, schemeDataOut ] = maskAndKeepMin(~, yIn, schemeDataIn)
% maskAndKeepMin: Example postTimestep processing routine.
%
%  [ yOut, schemeDataOut ] = maskAndKeepMin(t, yIn, schemeDataIn)
%
%  This function demonstrates two entirely different processes that
%    can be accomplished through the postTimestep odeCFLn integrator option.
%
%  The first is to mask the evolving implicit surface function
%    (or otherwise modify its value after each timestep).
%
%  In this case masking is accomplished by
%
%          yOut = max(yIn, schemeDataIn.mask);
%
%
%   which ensures that phi cannot be negative anywhere that mask is positive.
%
%  The second is to keep track of some feature of that implicit surface
%    function or otherwise modify the schemeData structure after each
%    timestep.
%
%  In this case the feature recorded is the pointwise minimum over time of phi
%
%          schemeDataOut.min = min(yIn, schemeDataIn.min);
%
%
% Parameters:
%   t              Current time.
%   yIn            Input version of the level set function, in vector form.
%   schemeDataIn   Input version of a structure (see below).
%
%   yOut           Output version of the level set function, in vector form.
%   schemeDataOut  Output version of the structure (possibly modified).
%
% schemeData is a structure containing data specific to this type of 
%   term approximation.  For this function it contains the field(s)
%
%   .doMask      Boolean specifying whether masking should be performed.
%   .doMin       Boolean specifying whether min should be taken.
%   .mask	 Function against which to mask the level set function.
%   .min         Function which stores the minimum of the level set
%                  function over time (it is modified at each timestep).
%
% schemeData may contain other fields.

  checkStructureFields(schemeDataIn, 'doMask', 'doMin');

  % Mask the current level set function.
  if(schemeDataIn.doMask)
    checkStructureFields(schemeDataIn, 'mask');
    yOut = max(yIn, schemeDataIn.mask);
  else
    yOut = yIn;
  end

  % Record any new minimum values for each node.
  %   Use yOut to get the masked version of the data (if masking).
  schemeDataOut = schemeDataIn;
  if(schemeDataIn.doMin)
    checkStructureFields(schemeDataIn, 'min');
    schemeDataOut.min = min(yOut, schemeDataOut.min);
  end
