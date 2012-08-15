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
%                  two dimensions.  Permits vector values for level parameter.
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
%    'wireframe' Show an isosurface as a wireframe.  Essentially the same
%                  as surface, but with the faces turned off and the edges
%                  turned on.

% Copyright 2007 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 6/29/04
% Modified to make use of gridnd2mesh, Ian Mitchell 5/17/07

%---------------------------------------------------------------------------
  if(nargin < 4)
    level = 0;
  end

  if((strcmp(displayType, 'contour') || strcmp(displayType, 'contourslice'))...
     && (prod(size(level)) == 1))
    % Scalar input to contour plot should be repeated.
    level = [ level level ];
  end

  if(~isempty(level))
    if((all(data(:) < min(level(:)))) | (all(data(:) > max(level(:)))))
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
    % In 2D, the visualization routines seem to be happy to use ndgrid.
    switch(displayType)
     case 'contour'
      [ garbage, h ] = contour(g.xs{1}, g.xs{2}, data, level, 'b');
      axis square;  axis manual;
     case 'surf'
      h = surf(g.xs{1}, g.xs{2}, data);
     otherwise
      error('Unknown display type %s for %d dimensional system', ...
            displayType, g.dim);    
    end
    
   %------------------------------------------------------------------------
   case 3
    % Stupid Matlab's stupid meshgrid vs ndgrid incompatibility really shows up
    % in 3D -- many of the 3D visualization routines only work for meshgrid
    % produced grids.  Therefore, we need to massage the grid and data to
    % make it work.
    [ mesh_xs, mesh_data ] = gridnd2mesh(g, data);
    
    switch(displayType)      
     case 'surface'
      h = patch(isosurface(mesh_xs{:}, mesh_data, level));
      isonormals(mesh_xs{:}, mesh_data, h);

      % This next line works without the conversion to meshgrid
      % coordinates; it appears to be the only 3D visualization code that
      % works without this conversion.
      %h = patch(isosurface(g.xs{1}, g.xs{2}, g.xs{3}, data, level));

      set(h, 'FaceColor', 'red', 'EdgeColor', 'none');

      % It would be nice to turn on a light here -- without it, the
      % surface always looks flat.  However, repeated calls to camlight
      % in the same axis produce repeated lights, which will quickly wash
      % out the pretty specular effects.  I can't find any easy way to
      % determine whether a light already exists for an axis.
      %camlight left;  camlight right;

      lighting phong;
      view(3)
     
     case 'slice'
      % For lack of a better idea of where to put the slices, we'll just slice
      % through the middle of the grid.  For final visualizations, users
      % should write their own code with a better choice of slice planes.
      avgx = mean(g.vs{1});
      avgy = mean(g.vs{2});
      avgz = mean(g.vs{3});
      h = slice(mesh_xs{:}, mesh_data, avgx, avgy, avgz);
      
     case 'contourslice'
      avgx = mean(g.vs{1});
      avgy = mean(g.vs{2});
      avgz = mean(g.vs{3});
      h = contourslice(mesh_xs{:}, mesh_data, ...
                       [ avgx ], [ avgy ], [ avgz ], level);
      %set(h, 'EdgeColor', 'black');
     
     case 'wireframe'
      % Because isosurface works on ndgrid, we don't need to use the converted
      % grid and data.
      h = patch(isosurface(g.xs{:}, data, level));
      set(h, 'FaceColor', 'none', 'EdgeColor', 'black');
      view(3)
 
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
