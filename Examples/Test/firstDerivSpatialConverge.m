% Script file to demonstrate the convergence rate of the first derivative
% approximation schemes.
%
% Uses firstDerivSpatialTest1 for a sequence of approximation schemes and
% grid sizes.
%
% Also demonstrates Matlab's bizarre but useful method of accessing
% structure fields by string variables.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
%
% This software is used, copied and distributed under the licensing 
% agreement contained in the file LICENSE in the top directory of 
% the distribution.
%
% Ian Mitchell, 1/25/04

run('../addPathToKernel');

% Be careful with grid dx that you don't overwhelm the memory
schemes = { @upwindFirstFirst; @upwindFirstENO2; ...
            @upwindFirstENO3a; @upwindFirstENO3b; ...
            @upwindFirstWENO5a; @upwindFirstWENO5b };
gridSizes = [ 20; 40; 80; 160; 320; 640; 1280 ];
%gridSizes = [ 20; 40; 80; 160 ];
dim = 2;
whichDim = 1;

% File id for a place to report the statistics (the screen).
screen = 1;

schemeNames = cell(size(schemes));
schemeStyles = { 'b-v', 'r-*', 'k-x', 'b-.o', 'k-+', 'r-.s' };

% Create some figures to display the results.
errorTypes = { 'maximum', 'average', 'rms', 'jumps' };
for errorTypeN = 1 : length(errorTypes)
  errorType = errorTypes{errorTypeN};
  f.(errorType) = figure;
  title([ 'Error analysis for ' errorType ' error.' ]);
  hold on;
end

h = zeros(length(errorTypes), length(schemes));
times = zeros(length(schemes), 1);
fprintf(screen, 'For grid size %d in %d dimensions, execution time (sec)\n',...
        gridSizes(end), dim);
for schemeN = 1 : length(schemes);
  scheme = schemes{schemeN};
  schemeNames{schemeN} = func2str(scheme);

  % Collect the error statistics.
  for gridSizeN = 1 : length(gridSizes);
    gridSize = gridSizes(gridSizeN);
    
    % We'll just look at the leftward error in the upwinded case.  It would
    % be nice to preallocate the errorStats array, but it is painfully
    % difficult to preallocate an array of structures.
    [ errorStats(gridSizeN), ignored, times(schemeN) ] ...
        = firstDerivSpatialTest1(scheme, dim, whichDim, 1.0 / gridSize); %#ok<SAGROW>
  end
  
  % Plot the error statistics.
  for errorTypeN = 1 : length(errorTypes)
    errorType = errorTypes{errorTypeN};
    figure(f.(errorType));
    h(errorTypeN, schemeN) = ...
        plot(gridSizes, [ errorStats(:).(errorType) ], schemeStyles{schemeN});
  end
  
  fprintf(screen, '\t%s\t%f\n', schemeNames{schemeN}, times(schemeN));
end

% Make the plots a little prettier.
for errorTypeN = 1 : length(errorTypes)
  errorType = errorTypes{errorTypeN};
  figure(f.(errorType));
  
  % Clean up the axes.
  set(gca, 'XScale', 'log', 'YScale', 'log', 'XTick', gridSizes, ...
           'XGrid', 'on', 'XMinorGrid', 'off', ...
           'YGrid', 'on', 'YMinorGrid', 'off');
  xlabel('Grid Size');  ylabel('Error');

  % Set the axes so that they include the data (plus some space vertically).
  currAxis = axis;
  newAxis = [ min(gridSizes), max(gridSizes), ...
              10^(floor(log10(currAxis(3)))), 10^(ceil(log10(currAxis(4)))) ];
  axis(newAxis);

  % Show slopes with various orders of convergence.
  slopePoint = [ gridSizes(4); 100 * newAxis(3) ];
  slopeHandles = addSlopes(slopePoint, gridSizes(5) - gridSizes(4), 'k-', ...
                           [ 0, -1, -2 -3 -5 ], ...
                           { '', '1', '2', '3', '5' });

  % Label the schemes.
  legend(h(errorTypeN, :), schemeNames, 3);

  % Make the lines and markers a bit more visible.
  set(h(errorTypeN, :), 'LineWidth', 2.0, 'MarkerSize', 9.0);
end
