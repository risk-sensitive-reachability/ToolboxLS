function [ data, g, data0 ] = tripleSine(b, accuracy)
% tripleSine: recreate figures 2.6 and 2.7 from Sethian.
%
%   [ data, g, data0 ] = tripleSine(b, accuracy)
%
% Recreates figures 2.6 and 2.7 from Sethian chapter 2, 
%   showing a combination of motion by surface normal and mean curvature
%   on a 2D example with a sinusoidal curve as the initial conditions.
%
% Basically, this is 2D motion in the normal direction with speed
%   			F(x) = 1 - b \kappa(x).
%
% The user may choose what multiplier to place on the curvature term
%   with input parameter b.
%
% The grid parameters have been chosen to try and get close to figure 2.6
%   with the same multiplier for \kappa; however, since no scaling
%   is given in the figures an exact match is not possible.
%  
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters:
%
%   b            Multiplier for the curvature dependent term.
%                  Must be nonnegative (default = 0).
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 4/21/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Default values.
if(nargin < 1)
  b = 0;
end

if(nargin < 2)
  accuracy = 'medium';
end

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 0.5;                  % End time.
plotSteps = 12;              % How many intermediate plots to produce?
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
useSubplots = 0;

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
g.min = [ 0; 0 ];
g.dx = 1 / 100;
g.max = [ +2; +2 ] - g.dx;
g.bdry = { @addGhostPeriodic; @addGhostExtrapolate };
g = processGrid(g);

%---------------------------------------------------------------------------
% What kind of display?
displayType = 'contour';

%---------------------------------------------------------------------------
% Create initial conditions (a sinusoidal wave).
%   Offset is designed to match Sethian's figures.
data = sin(3 * pi * (g.xs{1} + 0.5)) + 2 * g.xs{2} - 1.5;
data0 = data;

%---------------------------------------------------------------------------
% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
%   Same accuracy is used by both components of motion.
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

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Set up basic motion in the normal direction.
normalFunc = @termNormal;
normalData.grid = g;
normalData.speed = 1;
normalData.derivFunc = derivFunc;

%---------------------------------------------------------------------------
% Set up curvature motion.
curvatureFunc = @termCurvature;
curvatureData.grid = g;
curvatureData.curvatureFunc = @curvatureSecond;
curvatureData.b = b;

%---------------------------------------------------------------------------
% Combine components of motion.
if(b > 0)
  % If there is a nonzero curvature contribution to speed.
  schemeFunc = @termSum;
  schemeData.innerFunc = { normalFunc; curvatureFunc };
  schemeData.innerData = { normalData; curvatureData };
else
  % Otherwise ignore curvature.
  schemeFunc = normalFunc;
  schemeData = normalData;
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

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);
