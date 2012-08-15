function maxA = cellMatrixMax(A, takeAbs)
% cellMatrixMax: elementwise maximum over the state space of a matrix A(x)
%
%   maxA = cellMatrixMax(A, takeAbs)
%
% For a spatially varying matrix A stored as a cell matrix
%   (ie A{i,j} is an array which gives the value of A(i,j) at each node)
%   this function returns a regular matrix whose entries are the maximum
%   value over the spatial array (maximized independently in each element).
%
% parameters:
%   A		spatially varying matrix stored as a cell matrix
%   takeAbs	boolean, return max(abs(A{i,j})) instead of max(A{i,j})
%		  (optional, defaults to 0)
%
%   maxA	regular matrix containing the maximum entries
%
% so A(i,j) = max(A{i,j}) or max(abs(A{i,j}))

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 10/14/03

if(ndims(A) ~= 2)
  error('cellMatrixMax only works on 2D cell arrays');
end

if(nargin < 2)
  takeAbs = 0;
end

if(iscell(A))
  maxA = zeros(size(A));
  for i = 1 : size(A, 1)
    for j = 1 : size(A, 2)
      if(takeAbs)
        maxA(i,j) = max(abs(A{i,j}(:)));
      else
        maxA(i,j) = max(A{i,j}(:));
      end
    end
  end
else
  if(takeAbs)
    maxA = abs(A);
  else
    maxA = A;
  end
end
