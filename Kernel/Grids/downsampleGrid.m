function [ grid, varargout ] = downsampleGrid(grid, downsample_factor, varargin)
% downsampleGrid: Reduce the resolution of a grid.
%
%   [ grid, array1, array2, ... ] = downsampleGrid(grid, downsample_factor, array1, array2, ...)
%
% Reduces the number of grid nodes in each dimension by an integer amount,
% and suitably adjusts the upper bounds of the grid.  The lower bounds are
% not modified.  The resulting grid will have its number of nodes reduced by
% downsample_factor ^ grid.dim.
%
% Input Parameters:
%
%   grid: The grid structure to be downsampled.
%
%   downsample_factor: Positive integer.  The amount by which to
%   downsample.  Optional.  Default = 2.
%
%   array1, array2, ...: Arrays defined on the grid.  These arrays will
%   also be downsampled.  Optional.
%
% Output Parameters:
%
%   grid: The downsampled grid.
%
%   array1, array2, ...: The downsampled arrays.

% Copyright 2011 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Created by Ian Mitchell, 2011/01/30.
% Subversion tags for version control purposes.
% $Date: 2012-07-04 14:22:45 -0700 (Wed, 04 Jul 2012) $
% $Id: downsampleGrid.m 73 2012-07-04 21:22:45Z mitchell $

  %---------------------------------------------------------------------------
  %% Deal with optional but necessary input parameters.
  if(nargin < 2)
    downsample_factor = 2;
  end
  
  %---------------------------------------------------------------------------
  %% Index array for downsampling.
  index = cell(grid.dim, 1);
  for d = 1 : grid.dim
    index{d} = 1 : downsample_factor : grid.N(d);
  end
  
  %---------------------------------------------------------------------------
  %% Downsample any array arguments.
  % Do this downsampling first, since we need the old grid structure.
  for i = 1 : length(varargin)
    varargout{i} = varargin{i}(index{:});
  end
        
  %---------------------------------------------------------------------------
  %% Reduce the node position data.

  for d = 1 : grid.dim
    grid.vs{d} = grid.vs{d}(1 : downsample_factor : end);
    grid.xs{d} = grid.xs{d}(index{:});
  end

  %---------------------------------------------------------------------------
  %% Reduce the other fields as necessary.

  grid.dx = grid.dx * downsample_factor;
  
  % Need to be careful with N and max fields, since the downsample factor
  % may or may not include the previous last node.
  for d = 1 : grid.dim
    grid.N(d) = length(grid.vs{d});
    grid.max(d) = grid.vs{d}(end);
  end
  
  % Axis field is only defined for 2D and 3D grids; otherwise, it is empty.
  if ~isempty(grid.axis)
    grid.axis(2:2:end) = grid.max';
  end

  % Shape field is pretty standard except for one dimensional grids.
  if(grid.dim == 1)
    grid.shape = [ grid.N, 1 ];
  else
    grid.shape = grid.N';
  end

end

