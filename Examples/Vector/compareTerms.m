function [ dataC, dataH, g, data0 ] = ...
                       compareTerms(flowType, initShape, accuracy, displayType)
% compareTerms: compare convective and general Hamiltonian approximations
%
% [ dataC, dataH, g, data0 ] = ...
%                      compareTerms(flowType, initShape, accuracy, displayType)
%  
% This function demonstrates use of the vector level set methods by
%   comparing a upwinded convection approximation and a Lax-Friedrichs
%   general Hamiltonian approximation of rigid body motion according
%   to a constant vector field.
%
% The vector level set code is used in order to keep the two implicit surface
%   functions in temporal lockstep.  It is not really needed, but this
%   example provides a test of the code and an interesting comparison.
%
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters (all inputs have defaults):
%
%   flowType     String to specify type of flow field.
%                  'constant'    Constant flow field xdot = k (default).
%                  'linear'      Linear flow field xdot = A x.
%                  'converging'  A flow field which generates a shock.
%   initShape    String to specify what initial surface to use.
%	   	   'sphere'      A sphere in the appropriate dimension.
%	 	   'cube'  	 A cube in the appropriate dimension.
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
%   dataC        Solution at tMax using upwinded convection approximation.
%   dataH        Solution at tMax using Lax-Friedrichs Hamiltonian approx.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2005 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 8/30/04
% Modified for the version 1.1 vector level set formulation 3/03/05

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

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

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

% How many different term approximations do we have?
numElements = 2;

% Constants for each element.
elementC = 1;
elementH = 2;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 0;

%---------------------------------------------------------------------------
% Use periodic boundary conditions?
periodic = 0;

% Create the grid.
g.dim = 2;
g.min = -ones(g.dim, 1);
g.dx = 1 / 50;
if(periodic)
  g.max = (1 - g.dx) * ones(g.dim, 1);
  g.bdry = @addGhostPeriodic;
else
  g.max = ones(g.dim, 1);
  g.bdry = @addGhostExtrapolate;
end
g = processGrid(g);

%---------------------------------------------------------------------------
% Most of the time in constant flow case, we want flow in a 
%   distinguished direction, so assign first dimension's flow separately.
constantV = 0 * ones(g.dim);
constantV(1) = 2;
constantV = num2cell(constantV);

% Create linear flow field xdot = A * x
linearA = 2 * pi * [ 0 1 0 0; -1 0 0 0; 0 0 0 0; 0 0 0 0 ];      % Rotation.
%linearA = eye(4);                                               % Growth.
indices = { 1:g.dim; 1:g.dim };
linearV = cellMatrixMultiply(num2cell(linearA(indices{:})), g.xs);

% A flow field that creates a shock.
%   Constant flow in the first dimension.
%   All other dimensions flow towards zero.
shockV = cell(g.dim);
shockV{1} = +1;
for i = 2 : g.dim;
  shockV{i} = -sign(g.xs{i});
end

%---------------------------------------------------------------------------
if(nargin < 1)
  flowType = 'constant';
end

% Choose the flow field.
switch(flowType)
  
 case 'constant'
  v = constantV;
  
 case 'linear'
  v = linearV;

 case 'converging'
  v = shockV;
  
 otherwise
  error('Unknown flowType %s', flowType);
  
end

%---------------------------------------------------------------------------
% What kind of display?
if(nargin < 4)
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

% We will also need to differentiate the two implicit surface function plots.
switch(g.dim)
  case { 1, 2 }
    styleName = { 'LineStyle', 'LineWidth', 'Color' };
    styleC = { '-', 2, 'blue' };
    styleH = { '--', 2, 'red' };
  case 3
    styleName = { 'FaceColor', 'FaceAlpha' };
    styleC = { 'red', 0.6 };
    styleH = { 'blue', 0.6 };
  otherwise
    error('Display styles undefined for dimension %d', g.dim);
end

%---------------------------------------------------------------------------
if(nargin < 3)
  accuracy = 'low';
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
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Create initial conditions (a circle/sphere or square/cube).
%   Note that in the periodic BC case, these initial conditions will not
%   be continuous across the boundary unless the sphere/cube is centered.
%   In practice, we'll just ignore that little detail.
if(nargin < 2)
  initShape = 'sphere';
end

% These parameters can be used for either sphere or cube.
center = [ -0.4; 0.0; 0.0; 0.0 ];
radius = 0.4;

switch(initShape)
 case 'sphere'
  data0 = shapeSphere(g, center, radius);

 case 'cube'
  data0 = shapeRectangleByCenter(g, center, 2 * radius);

 otherwise
  error('Unknown initial shape string %s', initShape);
end

% Copy the initial conditions into a cell vector of implicit surface functions.
[ data{1 : numElements} ] = deal(data0);

%---------------------------------------------------------------------------
% Set up spatial approximation scheme for upwinded convection approximation.
schemeFunc{elementC} = @termConvection;
schemeData{elementC}.velocity = v;
schemeData{elementC}.grid = g;
schemeData{elementC}.derivFunc = derivFunc;

%---------------------------------------------------------------------------
% Set up spatial approximation scheme for Lax-Friedrichs approximation.
schemeFunc{elementH} = @termLaxFriedrichs;
schemeData{elementH}.grid = g;
schemeData{elementH}.hamFunc = @laxFriedrichsDemoHamFunc;
schemeData{elementH}.partialFunc = @laxFriedrichsDemoPartialFunc;
schemeData{elementH}.derivFunc = derivFunc;

% This parameter is not required by the termLaxFriedrichs code,
%   but is used by hamFunc and partialFunc.
schemeData{elementH}.velocity = v;                  

% Since the motion is convection, the type of artificial dissipation
%   is irrelevant.  Choose global since it is cheapest to compute.
dissType = 'global';

switch(dissType)
 case 'global'
  schemeData{elementH}.dissFunc = @artificialDissipationGLF;
 case 'local'
  schemeData{elementH}.dissFunc = @artificialDissipationLLF;
 case 'locallocal'
  schemeData{elementH}.dissFunc = @artificialDissipationLLLF;
 otherwise
  error('Unknown dissipation function %s', dissFunc);
end
  
%---------------------------------------------------------------------------
% Initialize Display
f = figure;
h = cell(numElements, 1);

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h{elementC} = visualizeLevelSet(g, data{elementC}, displayType, ...
                                level, [ 't = ' num2str(t0) ]);
setProperty(h{elementC}, styleName, styleC);
hold on;
h{elementH} = visualizeLevelSet(g, data{elementH}, displayType, ...
                                level, [ 't = ' num2str(t0) ]);
setProperty(h{elementH}, styleName, styleH);

if(g.dim > 1)
  axis(g.axis);
  daspect([ 1 1 1 ]);
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
y0 = cell(numElements, 1);
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  for i = 1 : numElements
    y0{i} = data{i}(:);
  end

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t, y, schemeData ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                               integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  for i = 1 : numElements
    data{i} = reshape(y{i}, g.shape);
  end

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
    hold on;
    if(g.dim > 1)
      axis(g.axis);
      daspect([ 1 1 1 ]);
    end
  end

  % Create new visualization.
  h{elementC} = visualizeLevelSet(g, data{elementC}, displayType, ...
                                  level, [ 't = ' num2str(tNow) ]);
  setProperty(h{elementC}, styleName, styleC);
  hold on;
  h{elementH} = visualizeLevelSet(g, data{elementH}, displayType, ...
                                  level, [ 't = ' num2str(tNow) ]);
  setProperty(h{elementH}, styleName, styleH);
  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);

dataC = data{elementC};
dataH = data{elementH};



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function hamValue = laxFriedrichsDemoHamFunc(~, data, deriv, schemeData)
% laxFriedrichsDemoHamFunc: demonstration analytic hamiltonian function.
%
% hamValue = laxFriedrichsDemoHamFunc(t, data, deriv, schemeData)
%
% This function implements the hamFunc prototype to demonstrate how
%   the Lax-Friedrichs HJ term approximation could be used to approximate
%   a convective flow field.  It calculates the analytic Hamiltonian
%   for such a flow field.
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
%   .velocity    Cell array containing the constant flow field.
%
% Ian Mitchell 2/11/04

checkStructureFields(schemeData, 'grid', 'velocity');

% Hamiltonian is just the dot product of velocity and costate.
hamValue = zeros(size(data));
for i = 1 : schemeData.grid.dim;
  hamValue = hamValue + schemeData.velocity{i} .* deriv{i};
end



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function alpha = ...
     laxFriedrichsDemoPartialFunc(~, ~, ~, ~, schemeData, dim)
% laxFriedrichsDemoPartialFunc: demonstration Hamiltonian partial function.
%
% alpha = ...
%    laxFriedrichsDemoPartialFunc(t, data, derivMin, derivMax, schemeData, dim)
%
% This function implements the partialFunc prototype to demonstrate how
%   the Lax-Friedrichs HJ term approximation could be used to approximate
%   a convective flow field.  It calculates the extrema of the absolute
%   value of the partials of the analytic Hamiltonian with respect to
%   the costate (gradient).
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
%   .velocity    Cell array containing the constant flow field.
%
% Ian Mitchell 2/11/04

checkStructureFields(schemeData, 'velocity');

% Constant flow field implies no dependence on costate.
alpha = abs(schemeData.velocity{dim});



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function setProperty(handle, propertyNames, propertyValues)
% setProperty: set handle graphics property settings.
%
%   setProperty(handle, propertyNames, propertyValues)
%
% For a given handle graphics object, sets one or more of that object's
%   graphics properties to particular values.
%
% Parameters:
%   handle          Graphics object handle(s).
%   propertyNames   Cell vector of names of the properties to be set.
%   propretyValues  Cell vector of values of the corresponding properties.

for i = 1 : length(propertyNames)
  for j = 1 : length(handle)
    set(handle(j), propertyNames{i}, propertyValues{i});
  end
end
drawnow;

