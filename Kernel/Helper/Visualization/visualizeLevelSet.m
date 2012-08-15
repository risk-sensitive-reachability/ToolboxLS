function h = visualizeLevelSet(g, data, displayType, level, titleStr)
% visualizeLevelSet: Display the level set at a particular time.
%
%   h = visualizeLevelSet(g, data, displayType, level, titleStr)
%
% Displays a variety of level set visualizations in dimensions 1 to 3.
%   The current figure and axis is used.
%
% A warning will be generated if the requested level set is missing.
%   For those display types that do not plot a level set
%   (such as a surf plot in 2D), the warning can be disabled by 
%   setting parameter level to be the empty vector [].
%
% Parameters:
%   g   	 Grid structure.
%   data         Array storing the implicit surface function.
%   displayType  String specifying the type of visualization (see below).
%   level        Which isosurface to display.  Defaults to 0.
%   titleStr     Optional string to place in the figure title.
%
%   h            Handle to the graphics object created.
%
% Display type options:
%
% Dimension 1:
%   'plot'       Plot the function value vs the state as a line.  If the 
%                  number of grid nodes is small, mark the node value as well.
%
% Dimension 2:
%   'contour'    Show an isocontour of the function as a solid curve in
%                  two dimensions.
%   'surf'       Plot the function value vs the two dimensional state as
%                  a surface plot.
%
% Dimension 3:
%    'surface'   Show an isosurface of the function as a solid surface in
%                  three dimensions.
%    'slice'     On slices through the x,y,z midpoints of the grid show
%                  the function value through color coding.
%    'contourslice'  On slices through the x,y,z midpoints of the grid
%                  show an isocontour of the function as a solid curve.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 6/29/04

%---------------------------------------------------------------------------
  if(nargin < 4)
    level = 0;
  end

  l = [ level level ];

  if(~isempty(level))
    if((all(data(:) < level)) | (all(data(:) > level)))
      warning('No implicitly defined surface exists');
    end
  end
  
%---------------------------------------------------------------------------
  switch(g.dim)
 
   %------------------------------------------------------------------------
   case 1
    switch(displayType)
     case 'plot'
      if(g.N < 20)
        % For very coarse grids, we can identify the individual nodes.
        h = plot(g.xs{1}, data, 'b-+');
      else
        h = plot(g.xs{1}, data, 'b-');
      end
     otherwise
      error('Unknown display type %s for %d dimensional system', ...
            displayType, g.dim);
    end
    
   %------------------------------------------------------------------------
   case 2
    switch(displayType)
     case 'contour'
      [ garbage, h ] = contour(g.xs{1}, g.xs{2}, data, l, 'b');
      axis square;  axis manual;
     case 'surf'
      h = surf(g.xs{1}, g.xs{2}, data);
     otherwise
      error('Unknown display type %s for %d dimensional system', ...
            displayType, g.dim);    
    end
    
   %------------------------------------------------------------------------
   case 3
    switch(displayType)
     case 'surface'
      h = patch(isosurface(g.xs{1}, g.xs{2}, g.xs{3}, data, level));
      %isonormals(g.xs{1}, g.xs{2}, g.xs{3}, data, h);
      set(h, 'FaceColor', 'red', 'EdgeColor', 'none');
      %camlight left;  camlight right;
      lighting phong;
      view(3)
     
     case 'slice'
      avgx = mean(g.vs{1});
      avgy = mean(g.vs{2});
      avgz = mean(g.vs{3});
      slice(g.xs{1}, g.xs{2}, g.xs{3}, data, avgx, avgy, avgz);
      
     case 'contourslice'
      avgx = mean(g.vs{1});
      avgy = mean(g.vs{2});
      avgz = mean(g.vs{3});
      h = contourslice(gx, gy, gz, ls, [ avgx ], [ avgy ], [ avgz ], l);
      %set(h, 'EdgeColor', 'black');
     
     otherwise
      error('Unknown display type %s for %d dimensional system', ...
            displayType, g.dim);
    end
    
   %------------------------------------------------------------------------
   otherwise
    warning('Unable to display data in dimension %d', g.dim);
    
  end
  
%---------------------------------------------------------------------------
  if(nargin >= 5)
    title(titleStr);
  end
  
  grid on;
  drawnow;
