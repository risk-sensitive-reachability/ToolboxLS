function [ errorL, errorR, time ] = ...
    firstDerivSpatialTest1(scheme, dim, whichDim, dx)
% firstDerivSpatialTest1: test various approximations of first derivative.
%
% [ errorL, errorR, time ] = firstDerivSpatialTest1(scheme, dim, whichDim, dx)
%
% Function to test the various approximations of the first spatial
% derivative.
%
% Uses a concatenation of shifted sin functions for data so as to introduce
% three derivative discontinuities (including one at the boundary in the
% periodic BC case).
%
% If the function passed in scheme begins with the characters "upwind"
% then both left and right approximations are tested, and both
% errorL and errorR are generated.
%
% Otherwise errorR = [] and errorL contains the error analysis.
%
% If no output parameters are requested, this function produces a figure
% showing the function and derivatives, as well as text output.
%
% Input Parameters:
%
%   scheme: Function handle to the derivative calculation scheme.
%
%   dim: How many dimensions should the computational grid be?
%
%   whichDim: Which dimension should the derivative be calculated in?
%
%   dx: Grid spacing.
%
% Output Parameters:
%
%   errorL: structure containing information about the error (error for the
%   left approximation in upwinded case).
%
%   errorR: structure containing information about the error for the right
%   approximation in the upwind case (otherwise just empty).
%
%   time: execution time to compute the derivative approximation
%
% The outputs errorL and errorR (when defined) are structures with the
% following fields:
%
%   .maximum    maximum error (inf norm)
%   .average    average error (one norm)
%   .rms        root mean square error (two norm)
%   .jump       average error at the three jumps

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 1/22/04

%---------------------------------------------------------------------------
% Some options.

% periodic BC?  Note that the analytic derivative calculation assumes periodic.
periodic = true;

% Is the derivative approximation upwinded?
schemeName = func2str(scheme);
if(strncmp('upwind', func2str(scheme), 6))
  upwinded = true;
else
  upwinded = false;
end

% A reasonably small number.
small = sqrt(eps);

% File id for a place to report the statistics (the screen).
screen = 1;

%---------------------------------------------------------------------------
% Build a grid.
g.dim = dim;
g.min = zeros(g.dim, 1);
g.dx = dx;
if(periodic)
  g.max = (1 - g.dx) * ones(g.dim, 1);
  g.bdry = @addGhostPeriodic;
else
  g.max = ones(g.dim, 1);
  g.bdry = @addGhostExtrapolate;
  g.bdryData.towardZero = 1;
end
g = processGrid(g);

%---------------------------------------------------------------------------
% Build a piecewise differentiable function.
xs = g.xs{whichDim};
xs2pi = 2 * pi * xs;

segments = 3;
indicator = cell(segments, 1);
indicator{1} = (xs - 0.25 < -small);
indicator{2} = ((xs - 0.25 > -small) & (xs - 0.5 < -small));
indicator{3} = (xs - 0.5 > -small);

data = (indicator{1} .* sin(xs2pi + 0.5 * pi) ...
        + indicator{2} .* sin(xs2pi - 0.5 * pi) ...
        + indicator{3} .* (sin(xs2pi) + 1));

%---------------------------------------------------------------------------
% Compute the analytic derivatives.
if(~periodic)
  warning('Analytic derivative approximation assumes periodic BC.');
end

% A cell array of indices which shifts dimension of interest by 1 element.
shiftBy1 = cell(g.dim, 1);
for i = 1 : g.dim;
  shiftBy1{i} = 1 : g.N(i);
end
shiftBy1{whichDim} = [ g.N(whichDim), 1 : g.N(whichDim) - 1 ];

analyticL = 2 * pi * (indicator{1}(shiftBy1{:}) .* cos(xs2pi + 0.5 * pi) ...
                      + indicator{2}(shiftBy1{:}) .* cos(xs2pi - 0.5 * pi) ...
                      + indicator{3}(shiftBy1{:}) .* cos(xs2pi));
analyticR = 2 * pi * (indicator{1} .* cos(xs2pi + 0.5 * pi) ...
                      + indicator{2} .* cos(xs2pi - 0.5 * pi) ...
                      + indicator{3} .* cos(xs2pi));
analytic0 = 0.5 * (analyticL + analyticR);

%---------------------------------------------------------------------------
% What is the true derivative (at the jumps, keep both values)?
derivX = cell(segments, 1);
derivX{1} = linspace(0, 0.25, 100)';
derivX{2} = linspace(0.25, 0.5, 100)';
derivX{3} = linspace(0.5, 1, 200)';
derivY = cell(segments, 1);
derivY{1} = 2 * pi * cos(2 * pi * derivX{1} + 0.5 * pi);
derivY{2} = 2 * pi * cos(2 * pi * derivX{2} - 0.5 * pi);
derivY{3} = 2 * pi * cos(2 * pi * derivX{3});

%---------------------------------------------------------------------------
% Calculate the derivative approximation(s)

startTime = cputime;
if(upwinded)
  [ derivL, derivR ] = feval(scheme, g, data, whichDim);
else
  deriv0 = feval(scheme, g, data, whichDim);
end
endTime = cputime;

%---------------------------------------------------------------------------
if(nargout == 0)
  % Display the analytic results.
  figure;
  hold on;
  
  % We need to get just a 1D slice of the data set.
  indices = num2cell(ones(g.dim, 1));
  indices{whichDim} = 1 : g.N(whichDim);
  
  % Plot the function.
  h(1) = plot(squeeze(xs(indices{:})), squeeze(data(indices{:})), 'm:');
  
  % Plot the analytic derivative.
  for i = 1 : segments
    h(2) = plot(derivX{i}, derivY{i}, 'k--');
  end
  
  %---------------------------------------------------------------------------
  % Compute and display the derivative approximation.
  if(upwinded)
    h(3) = plot(squeeze(xs(indices{:})), squeeze(derivL(indices{:})), 'b-*');
    h(4) = plot(squeeze(xs(indices{:})), squeeze(derivR(indices{:})), 'r-*');
    legend(h, { 'function'; 'analytic D'; 'D-'; 'D+' });
  else
    h(3) = plot(squeeze(xs(indices{:})), squeeze(deriv0(indices{:})), 'b-*');
    legend(h, { 'function'; 'analytic D'; 'D0' });
  end
  
  title([ 'Comparing derivative approximations with ' func2str(scheme) ]);
  % This is why we don't use 'g' instead of 'grid' in the scripts
  grid on;
  
  %---------------------------------------------------------------------------
  % Report the statistics.
  fprintf(screen, 'Computational time %f\n', endTime - startTime);
  meanD = mean(abs(analytic0(:)));
  fprintf(screen, 'Error relative to mean magnitude derivative %g\n', meanD);
  if(upwinded)
    errorL = abs(analyticL - derivL);
    errorR = abs(analyticR - derivR);
    fprintf(screen, 'Relative left error  max %g, mean %g\n', ...
            max(errorL(:)) / meanD, mean(errorL(:)) / meanD);
    fprintf(screen, 'Relative right error max %g, mean %g\n', ...
            max(errorR(:)) / meanD, mean(errorR(:)) / meanD);
  else
    error0 = abs(analytic0 - deriv0);
    fprintf(screen, 'Relative centered Error  max %g, mean %g\n', ...
            max(error0(:)) / meanD, mean(error0(:)) / meanD)
  end

%---------------------------------------------------------------------------
else

  jumps = [ 0; 0.25; 0.5 ];
  
  % Create cell array with array indices so as to strip out a single
  %   column vector slice of the data array in the right dimension
  indices = num2cell(ones(g.dim, 1));
  indices{whichDim} = 1 : g.N(whichDim);
  xs = squeeze(g.xs{whichDim}(indices{:}));
  
  % Calculate statistics and return them.
  if(upwinded)
    errorVector = analyticL - derivL;
    errorVector = squeeze(errorVector(indices{:}));
    errorL = analyze(xs, errorVector, jumps, small);
    errorVector = analyticR - derivR;
    errorVector = squeeze(errorVector(indices{:}));
    errorR = analyze(xs, errorVector, jumps, small);    
    
  else
    errorVector = analytic0 - deriv0;
    errorVector = squeeze(errorVector(indices{:}));
    errorL = analyze(xs, errorVector, jumps, small);    
    errorR = [];
  end
  
  time = endTime - startTime;
end



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function errorStruct = analyze(xs, errorVector, jumps, small)
%  errorStruct = analyze(xs, errorVector, jumps, small)
%
% Helper function to collect statistics about the error.

n = length(errorVector);
  
errorStruct.maximum = norm(errorVector, inf);
errorStruct.average = norm(errorVector, 1) / n;
errorStruct.rms = norm(errorVector, 2) / sqrt(n);

errorJump = zeros(length(jumps), 1);
for i = 1 : length(jumps);

  jumpX = jumps(i);
  jumpI = find(abs(xs - jumpX) < small);
  errorJump(i) = errorVector(jumpI);
  
end

errorStruct.jumps = mean(abs(errorJump));
