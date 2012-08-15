function [ ydot, stepBound ] = termLaxFriedrichs(t, y, schemeData)
% termLaxFriedrichs: approximate H(x,p) term in an HJ PDE with Lax-Friedrichs.
%
% function [ ydot, stepBound ] = termLaxFriedrichs(t, y, schemeData)
%
% Computes a Lax-Friedrichs (LF) approximation of a general Hamilton-Jacobi
%   equation.  Global LF, Local LF, Local Local LF and Stencil LF are
%   implemented by choosing different dissipation functions.  The PDE is:
%
%            D_t \phi = -H(x, D_x \phi).
%
% Based on methods outlined in O&F, chapter 5.3 and 5.3.1.
%
% Parameters:
%   t            Time at beginning of timestep.
%   y            Data array in vector form.
%   schemeData	 A structure (see below).
%
%   ydot	 Change in the data array, in vector form.
%   stepBound	 CFL bound on timestep for stability.
%
% schemeData is a structure containing data specific to this type of 
%   term approximation.  For this function it contains the field(s)
%
%   .grid	 Grid structure (see processGrid.m for details).
%   .derivFunc   Function handle to upwinded finite difference 
%                  derivative approximation.
%   .dissFunc    Function handle to LF dissipation calculator.
%   .hamFunc	 Function handle to analytic hamiltonian H(x,p).
%   .partialFunc Function handle to extrema of \partial H(x,p) / \partial p.
%
% Note that options for derivFunc and dissFunc are provided as part of the
%   level set toolbox, while hamFunc and partialFunc depend on the exact
%   term H(x,p) and are user supplied.  Note also that schemeData may
%   contain addition fields at the user's discretion; for example, fields
%   containing parameters useful to hamFunc or partialFunc.
%
%
% schemeData.hamFunc should have prototype
%
%         hamValue = hamFunc(t, data, deriv, schemeData)
%
%   where t and schemeData are passed directly from this function, data = y
%   has been reshaped into its original size, and deriv is a cell vector (of
%   length grid.dim) containing the elements of the costate p = \grad \phi.
%   The return value should be an array (the size of data) containing
%   H(x,p).
%
%
% For details on schemeData.partialFunc, see the dissipation functions.
%
%
% In the notation of OF text:
%
%   data	  \phi.
%   derivFunc	  Function to calculate phi_i^+-.
%   dissFunc      Function to calculate the terms with alpha in them.
%   hamFunc	  Function to calculate analytic H.
%   partialFunc	  \alpha^i (dimension i is an argument to partialFunc).
%
%   update	  -\hat H.


% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell 5/13/03
% Calling parameters significantly modified, Ian Mitchell 2/11/04.

  %---------------------------------------------------------------------------
  checkStructureFields(schemeData, 'grid', 'derivFunc', 'dissFunc', ...
                       'hamFunc', 'partialFunc');
    
  %---------------------------------------------------------------------------
  grid = schemeData.grid;
  data = reshape(y, grid.shape);

  %---------------------------------------------------------------------------
  % Get upwinded and centered derivative approximations.
  derivL = cell(grid.dim, 1);
  derivR = cell(grid.dim, 1);
  derivC = cell(grid.dim, 1);
  
  for i = 1 : grid.dim
    [ derivL{i}, derivR{i} ] = feval(schemeData.derivFunc, grid, data, i);
    derivC{i} = 0.5 * (derivL{i} + derivR{i});
  end

  %---------------------------------------------------------------------------
  % Calculate update: (unstable) analytic hamiltonian
  %                   - (dissipative) stabiliziation.

  ham = feval(schemeData.hamFunc, t, data, derivC, schemeData);
  [ diss, stepBound ] = feval(schemeData.dissFunc, t, data, ...
                              derivL, derivR, schemeData);
  
  delta = ham - diss;
  
  %---------------------------------------------------------------------------
  % Reshape output into vector format and negate for RHS of ODE.
  ydot = -delta(:);
