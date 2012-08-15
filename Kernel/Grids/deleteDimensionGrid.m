function grid = deleteDimensionGrid(grid, dimensions)
% deleteDimensionGrid: Project away one or more dimensions in a grid.
%
%   grid = deleteDimensionGrid(grid, dimensions)
%
% Reduces the dimension of a grid by projecting away one or more of the
% dimensions.  Those dimensions are removed entirely, and the data in the
% grid.xs field is appropriately modified.
%
% Input Parameters:
%
%   grid: The grid structure to be modified.
%
%   dimensions: Positive integer or vector of positive integers.
%   Indexes of dimension(s) to be removed.
%
% Output Parameters:
%
%   grid: The grid without the requested dimension.

% Copyright 2011 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Created by Ian Mitchell, 2011/01/30.
% Subversion tags for version control purposes.
% $Date: 2011-05-16 16:06:25 -0700 (Mon, 16 May 2011) $
% $Id: deleteDimensionGrid.m 66 2011-05-16 23:06:25Z mitchell $

  %---------------------------------------------------------------------------
  %% What we actually need is a vector of the dimensions to keep.
  % I don't know of a way to vectorize this operation, but this loop will
  % be so cheap compared to the recreation of grid.xs that it hardly
  % matters.
  keep_dims = zeros(grid.dim - numel(dimensions), 1);
  count = 1;
  for i = 1 : grid.dim;
    if all(i ~= dimensions)
      keep_dims(count) = i;
      count = count + 1;
    end
  end
  
  %---------------------------------------------------------------------------
  %% The grid is lower dimensional.
  grid.dim = grid.dim - numel(dimensions);

  %---------------------------------------------------------------------------
  %% Remove the desired dimensions from all fields that are simple vectors.
  grid.min = grid.min(keep_dims);
  grid.max = grid.max(keep_dims);
  grid.N = grid.N(keep_dims);
  grid.dx = grid.dx(keep_dims);
  
  grid.vs = grid.vs(keep_dims);
  grid.bdry = grid.bdry(keep_dims);
  grid.bdryData = grid.bdryData(keep_dims);

  %---------------------------------------------------------------------------
  %% For remaining fields, it is easier to recreate them.
  
  grid.xs = cell(grid.dim, 1);
  if(grid.dim > 1)
    [ grid.xs{:} ] = ndgrid(grid.vs{:});
  else
    grid.xs{1} = grid.vs{1};
  end
  
  % Axis field is only defined for 2D and 3D grids; otherwise, it is empty.
  if((grid.dim == 2) || (grid.dim == 3))
    grid.axis = zeros(1, 2*grid.dim);
    grid.axis(1:2:end) = grid.min';
    grid.axis(2:2:end) = grid.max';
  end

  % Shape field is pretty standard except for one dimensional grids.
  if(grid.dim == 1)
    grid.shape = [ grid.N, 1 ];
  else
    grid.shape = grid.N';
  end

end % deleteDimensionGrid().
