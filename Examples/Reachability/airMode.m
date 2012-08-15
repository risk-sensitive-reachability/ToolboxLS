function [ reach, g, avoid, data0 ] = airMode(accuracy)
% airMode: demonstrate the 3 mode collision avoidance scenario.
%
%   [ reach, g, avoid, data0 ] = airMode(accuracy)
%  
% In this example, the target set is a circle at the origin
%   that represents a collision in relative coordinates between the evader
%   (player a, fixed at the origin facing right) and the pursuer (player b).
%
% Both aircraft are flying a fixed, synchronized series of segments.
%   Each segment has fixed linear and angular velocity.  The only input
%   is the time at which the protocol is initiated, which is a controlled
%   switch from mode 1 to mode 2.  The switch from mode 2 to mode 3 is timed.
%
% The relative coordinate dynamics in mode 1 and 3 are
%
%        \dot x    = -v_a + v_b \cos \psi
%	 \dot y    = v_b \sin \psi
%
% The relative coordinate dynamics in mode 2 are
%
%        \dot x    = -v_a + v_b \cos \psi + \omega y
%	 \dot y    = v_b \sin \psi - \omega x
%
%   where \omega, \psi, v_a and v_b are constants.
%
% For more details, see my PhD thesis, section 6.1.
%
% Because there are no continuous inputs, this reachable set calculation
%   can be performed with convection instead of general HJ PDEs.
%
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   aircraft parameters, etc.
%
% To get exactly the result from the thesis choose:
%   targetRadius = 5, velocityA =  3, velocityB = 4, 
%   \psi = 4 \pi / 3, \omega = 1.
%
% Parameters:
%
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%
%   reach        Implicit surface function for unsafe (reach) set in mode 1.
%   g            Grid structure on which data was computed.
%   avoid        Implicit surface function for safe to switch set (avoid)
%                  in mode 1.
%   data0        Implicit surface function for target set.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 4/18/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Problem Parameters.
%   targetRadius  Radius of target circle (positive).
%   velocityA	  Speed of the evader (positive constant).
%   velocityB	  Speed of the pursuer (positive constant).
%   psi           Relative heading of the two vehicles.
%   omega         Angular velocity in mode 2.
targetRadius = 5;
velocityA = 3;
velocityB = 4;
psi = -4 * pi / 3;
omega = 1;

%---------------------------------------------------------------------------
% Integration parameters.
tMaxStraight = 5;                    % End time for straight segments.
tMaxCurved = pi / omega;             % End time for curved segments.

% Period at which intermediate plots should be produced.
tPlot = 1;

% A figure to put intermediate plots.
fig = figure;

%---------------------------------------------------------------------------
% Default accuracy.
if(nargin < 1)
  accuracy = 'medium';
end

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
g.min = -25;
g.max = +25;
g.bdry = @addGhostExtrapolate;
g.N = 101;
g = processGrid(g);

%---------------------------------------------------------------------------
% Create initial conditions (circle centered on origin).
data0 = shapeSphere(g, [ 0; 0 ], targetRadius);

%---------------------------------------------------------------------------
% Create the flow fields for the two types of motion.
%   Multiply by -1 to get forward PDE.
straight = { -(-velocityA + velocityB * cos(psi)); ...
             -(velocityB * sin(psi)) };
curved = { -(-velocityA + velocityB * cos(psi) + omega * g.xs{2}); ...
           -(velocityB * sin(psi) - omega * g.xs{1}) };

%---------------------------------------------------------------------------
% Mode 3 reachable set is a simple reach computation.
%   Anywhere in this set in mode 3 is unsafe.
[ mode3, stepTime ] = findReachSet(g, data0, straight, accuracy, ...
                                   tMaxStraight, tPlot, fig, 1, []);
executionTime = stepTime;

%---------------------------------------------------------------------------
% Next is Mode 2 reachable set without relation to mode 3.
%   Anywhere in this set in mode 2 is unsafe.
[ mode2, stepTime ] = findReachSet(g, data0, curved, accuracy, ...
                                   tMaxCurved, tPlot, fig, 1, []);
executionTime = executionTime + stepTime;

%---------------------------------------------------------------------------
% Now we need the set of states that is taken by mode 2 into
%   an unsafe state in mode 3.
%   This set is time dependent, but we are only interested in its
%   shape at the time of the switch from mode 1.

% Initial condition for this set are taken from mode 3 
%   (taking rotational reset map between modes 2 and 3 into account).
initial = rot90(mode3, -1);

[ switch2, stepTime ] = findReachSet(g, initial, curved, accuracy, ...
                                     tMaxCurved, tPlot, fig, 0, []);
executionTime = executionTime + stepTime;

%---------------------------------------------------------------------------
% Finally, we get back to mode 1, which has a controlled switch and hence
%   an escape set.

% Escape set is complement of 
%   everything unsafe in mode 2 
%   union everything mode 2 will take to unsafety in mode 3.
%   (don't forget the reset map between modes 1 and 2.)
avoid = min(rot90(mode2, -1), rot90(switch2, -1));

[ reach, stepTime ] = findReachSet(g, data0, straight, accuracy, ...
                                   tMaxStraight, tPlot, fig, 1, avoid);
executionTime = executionTime + stepTime;

%---------------------------------------------------------------------------
% Now we can summarize the results in mode 1.
figure;
level = [ 0 0 ];
contourf(g.xs{1}, g.xs{2}, -reach, level);
hold on;
contour(g.xs{1}, g.xs{2}, avoid, level, 'r-');
contour(g.xs{1}, g.xs{2}, mode3, level, 'b--');
axis equal;
axis(g.axis);
fprintf('Total execution time %g seconds', executionTime);

% Clip back the axis bounds slightly so that the distortion caused
%   by the artificial boundary doesn't show.
clip = [ +5, -5, +5, -5 ];
axis(g.axis + clip);


%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function [ final, executionTime ] = findReachSet(g, initial, velocity, ...
                                accuracy, tMax, tPlot, fig, growOnly, avoid)
% findReachSet: Compute reach or reach-avoid sets for multimode protocol.
%
% [ final, executionTime ] = findReachSet(g, initial, velocity, ...
%                                accuracy, tMax, tPlot, fig, growOnly, avoid)
%
% Assembles the schemeFunc data structure and runs the HJ PDE to compute
%   the reach or reach-avoid set for a single mode of a multimode collision
%   avoidance protocol.  Specialized to handle convective dynamics
%   within each individual mode
%   (ie termConvection rather than termLaxFriedrichs).
%
% Parameters:
%
%   g             Grid structure on which data was computed.
%   initial       Array containing the implicit surface for the target set.
%   velocity      Cell vector containing the convective flow field.
%   accuracy      Controls the order of approximations.
%                   'low'         Use odeCFL1 and upwindFirstFirst.
%                   'medium'      Use odeCFL2 and upwindFirstENO2.
%                   'high'        Use odeCFL3 and upwindFirstENO3.
%                   'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   tMax          Integration time (integrations start at t0 = 0).
%   tPlot         Time at which to produce intermediate figures.
%   fig           Handle to figure in which to produce intermediate figures.
%   growOnly      Boolean, specifies whether to restrict the temporal
%                   derivative of the implicit surface function so that the
%                   reachable set only grows.
%   avoid         Array containing the implicit surface for the escape set.
%                   Set to [] if there is no escape set.
%
%
%   final         Array containing the implicit surface function of the
%                   reach or reach-avoid set.
%   executionTime Time required to compute final (from cputime, in seconds).

% Ian Mitchell, 4/18/04

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Visualize the 2D reachable set.
displayType = 'contour';

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 1;

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

t0 = 0;                              % Start time.

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.grid = g;
schemeData.velocity = velocity;

%---------------------------------------------------------------------------
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

%---------------------------------------------------------------------------
% Restrict the Hamiltonian so that reachable set only grows.

if(growOnly)
  innerFunc = schemeFunc;
  innerData = schemeData;
  clear schemeFunc schemeData;

  % Wrap the true Hamiltonian inside the term approximation restriction.
  schemeFunc = @termRestrictUpdate;
  schemeData.innerFunc = innerFunc;
  schemeData.innerData = innerData;
  schemeData.positive = 0;
end

%---------------------------------------------------------------------------
if(isempty(avoid))
  data = initial;
else
  % Ensure that the initial data satisfies the avoid set.
  data = max(initial, avoid);

  % Set up data required for masking by the avoid set.
  %   Mask will be compared to vector form of data array used by integrator.
  schemeData.maskData = avoid(:);
  schemeData.maskFunc = @max;

  % Let the integrator know what function to call.
  integratorOptions = odeCFLset(integratorOptions, ...
                                'postTimestep', @postTimestepMask);
end

%---------------------------------------------------------------------------
% Initialize Display
figure(fig);

h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);

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

  % Get correct figure.
  figure(fig);

  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);

end

executionTime = cputime - startTime;
fprintf('Execution time %g seconds\n', executionTime);

final = data;
