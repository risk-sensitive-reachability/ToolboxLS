function [ data, g, data0 ] = acoustic(filename, accuracy, compress)
% animateAcoustic: create an animation of the growth of the acoustic reach set.
%
%   [ data, g, data0 ] = animateAcoustic(filename, accuracy, compress)
%  
% This file generates an animation showing how the reachable set grows
%   as time progresses.  It is basically a combination of the files:
%        acoustic (which sets up and handles the reach set calculation),
%        spinAnimation (which has the code necessary for animations).
%
% NOTES:
%
% 1) This function will probably only work in Windows, since it uses
%    the avifile command and the avi animation format.
%
% 2) While the animation is being generated, don't move anything
%    (mouse, other windows) in front of the matlab figure window,
%    otherwise that image will be captured into the avi.
%
% 3) If this file stops because of an error, the avi file will remain open.
%    In that case, you must issue a "clear all" command to get rid of
%    the open file handle and then delete the partially completed avi file.
%
%  
% In this example the target set is a horizontal wide but shallow rectangle
%   near the origin, which represents the pursuer's capture set.
%   The computation is done in relative coordinates with the pursuer fixed at
%   the origin.
%
% The relative coordinate dynamics are
%
%        \dot x    = W_p ( 0) + W_p/R ( y) b + 2 W_e S_e ax
%	 \dot y    = W_p (-1) + W_p/R (-x) b + 2 W_e S_e ay
%
%   where W_p, W_e, R, S are constants
%         S_e = min(sqrt(x^2 + y^2), S)
%         input a is trying to avoid the target, sqrt(ax^2 + ay^2) \leq 1.
%	  input b is trying to hit the target, |b| \leq 1.
%
% For more details, see my PhD thesis, section 3.2.
%
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   aircraft parameters, etc.
%
% To get exactly the result from the thesis choose:
%   We = 1.3, Wp = 1, R = 0.8, S = 0.5
%   and capture set [ -3.5, +3.5 ] x [ -0.2, 0 ]
%
% Parameters:
%
%   filename     The name to give to the animation avi file.
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst.
%                  'medium'      Use odeCFL2 and upwindFirstENO2 (default).
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   compress     Boolean specifying whether to use lossy compression to
%                  (significantly) reduce the file size.  The degree of
%                  compression can be modified by changing the source code.
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 7/14/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Some animation parameters.
frames = 131;
deleteLabels = 0;               % Labels tend to float around and look bad.
qualityValue = 90;		% Even 100 gets significant compression.

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 3.9;                  % End time.
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (frames - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

% What kind of dissipation?
dissType = 'global';

%---------------------------------------------------------------------------
% Problem Parameters.
%   We         Speed of evader.
%   Wp         Speed of pursuer.
%   R          Turn radius of pursuer.
%   S          Radius of maximum speed for evader.
%   boundA     Range of norm of input A (evader).
%   boundB     Range of magnitude of input B (pursuer).
%   rangeTx    Range of target set in x direction.
%   rangeTy    Range of target set in y direction.
We = 1.3;
Wp = 1.0;
R = 0.8;
S = 0.5;
boundA = 1.0;
boundB = 1.0;
rangeTx = [ -3.5; 3.5 ];
rangeTy = [ -0.2; 0.0 ];

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Visualize the 3D reachable set.
displayType = 'contour';

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
g.min = [ -4; -2 ];
g.max = [ +4; +3 ];
g.bdry = @addGhostExtrapolate;
g.dx = 1/25;
g = processGrid(g);

%---------------------------------------------------------------------------
% Create initial conditions.
%   A rectangle is the intersection of four halfplanes.
data = shapeRectangleByCorners(g, [ rangeTx(1); rangeTy(1) ], ...
                                  [ rangeTx(2); rangeTy(2) ]);
data0 = data;

%---------------------------------------------------------------------------
% Set up spatial approximation scheme.
schemeFunc = @termLaxFriedrichs;
schemeData.hamFunc = @acousticHamFunc;
schemeData.partialFunc = @acousticPartialFunc;
schemeData.grid = g;

% The Hamiltonian and partial functions need problem parameters.
schemeData.We = We;
schemeData.Wp = Wp;
schemeData.R = R;
schemeData.boundA = boundA;
schemeData.boundB = boundB;

% For evaluating the evader's speed, 
%   we might as well precompute the term min(sqrt(x^2 + y^2), S).
schemeData.speedBound = min(sqrt(g.xs{1}.^2 + g.xs{2}.^2), S);

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
if(nargin < 1)
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
% Initialize Display
fig = figure;
h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);
axis(g.axis);

%---------------------------------------------------------------------------
% Now establish the animation stuff.

% Set up the figure window nicely (you can change the resolution here).
set(fig, 'Position', [ 100 100 384 384 ], 'color', 'white');
box on

% Turn off stuff we don't want -- animations are typically too busy looking.
if(deleteLabels)
  set(gca, 'XTick', [], 'YTick', []);
  delete(get(gca, 'Title'));
  delete(get(gca, 'XLabel'));
  delete(get(gca, 'YLabel'));
end

% Create the avi file (choose a smaller qualityValue to get smaller files).
if((nargin < 3) | compress)
  mov = avifile(filename, 'quality', qualityValue);
else
  mov = avifile(filename, 'compression', 'none');
end

%---------------------------------------------------------------------------
% Capture the first frame.

frame = 1;
drawnow;
f = getframe(gcf);
mov = addframe(mov, f);

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

  % Delete last visualization if necessary.
  delete(h);

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);
  axis(g.axis);

  % Get next frame.
  frame = frame + 1;
  drawnow;
  f = getframe(gcf);
  mov = addframe(mov, f);

end

endTime = cputime;
fprintf('Total execution time %g seconds', endTime - startTime);

% We're finished with the movie.
mov = close(mov);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = acousticHamFunc(t, data, deriv, schemeData)
% acousticHamFunc: analytic Hamiltonian for acoustic capture game.
%
% hamValue = acousticHamFunc(t, data, deriv, schemeData)
%
% This function implements the hamFunc prototype for the acoustic
%   capture game.  The Hamiltonian can be found in my PhD thesis
%   chapter 3.2, bottom of page 57.
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
%   .We          Speed of evader.
%   .Wp          Speed of pursuer.
%   .R           Turn radius of pursuer.
%   .boundA      Range of norm of input A (evader).
%   .boundB      Range of magnitude of input B (pursuer).
%   .speedBound  precomputed term min(sqrt(x^2 + y^2), S).
%
% Ian Mitchell 4/21/04

checkStructureFields(schemeData, 'grid', 'We', 'Wp', 'R', 'speedBound', ...
                                 'boundA', 'boundB');

grid = schemeData.grid;

% implements equation at the bottom of p.57 from my thesis term by term
%   with allowances for nonunit \script A and \script B
%   where deriv{i} is p_i
%         x is grid.xs{1}, y is grid.xs{2}
%         \script A is boundA and \script B is boundB
hamValue = -(-schemeData.Wp * deriv{2} ...
             -schemeData.Wp / schemeData.R * schemeData.boundB ...
              * abs(deriv{1} .* grid.xs{2} - deriv{2} .* grid.xs{1}) ...
             + 2 * schemeData.We * schemeData.boundA * schemeData.speedBound...
               .* sqrt(deriv{1}.^2 + deriv{2}.^2));



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = ...
              acousticPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
% acousticPartialFunc: Hamiltonian partial fcn for acoustic capture game.
%
% alpha = acousticPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
%
% This function implements the partialFunc prototype for the acoustic
%   capture game.  The Hamiltonian can be found in my PhD thesis
%   chapter 3.2, bottom of page 57.
%
% It calculates the extrema of the absolute value of the partials of the 
%   analytic Hamiltonian with respect to the costate (gradient).
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
%
%   .grid	 Grid structure.
%   .We          Speed of evader.
%   .Wp          Speed of pursuer.
%   .R           Turn radius of pursuer.
%   .boundA      Range of norm of input A (evader).
%   .boundB      Range of magnitude of input B (pursuer).
%   .speedBound  precomputed term min(sqrt(x^2 + y^2), S).
%
% Ian Mitchell 4/21/04

checkStructureFields(schemeData, 'grid', 'We', 'Wp', 'R', 'speedBound', ...
                                 'boundA', 'boundB');

grid = schemeData.grid;

% To bound the effect of the evader, we assume for each value of dim
%   that the evader's entire input is in that direction alone.
evaderSpeed = 2 * schemeData.We * schemeData.boundA * schemeData.speedBound;

switch dim
  case 1
    alpha = (schemeData.Wp / schemeData.R ...
             * schemeData.boundB * abs(grid.xs{2}) + evaderSpeed);

  case 2
    alpha = (schemeData.Wp + schemeData.Wp / schemeData.R ...
             * schemeData.boundB * abs(grid.xs{1}) + evaderSpeed);

  otherwise
    error([ 'Partials for the acoustic capture game' ...
            ' only exist in dimensions 1-2' ]);
end
