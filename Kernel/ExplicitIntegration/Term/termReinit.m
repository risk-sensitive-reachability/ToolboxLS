function [ ydot, stepBound ] = termReinit(t, y, schemeData)
% termReinit: a Godunov solver for the reinitialization HJ PDE.
%
% [ ydot, stepBound ] = termReinit(t, y, schemeData)
%
% Computes a Godunov approximation to motion by the reinitialization
%   equation.  While the reinitialization equation is a general nonlinear HJ
%   PDE, such a Godunov approximation is the least dissipative monotone
%   approximation (less dissipative than Roe-Fix or Lax-Friedrichs).  The
%   reinitialization equation is
%
%            D_t \phi = -sign(\phi_0)(\|\grad \phi\| - 1).
%
%   where phi_0 is the initial conditions.  Solving the reinitialization
%   equation turns an implicit surface function into a signed distance
%   function.  It is iterative, and often slower than a fast marching
%   method; however, it can use high order approximations and can start
%   directly from the implicit surface function.
%
% The reinitialization equation is discussed in O&F chapter 7.4.  The
%   Gudonov solver used below comes from from appendix A.3 of Fedkiw, Aslam,
%   Merriman & Osher, JCP 152, pp. 457-492 (1999) (citation [63] in O&F).
%
% The sign() approximation is given in the subfunction smearedSign.
%
% parameters:
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
%   .initial	 initial implicit surface function
%                (used to determine on which side of surface node should lie)
%
% It may contain addition fields at the user's discretion.
%
% In the notation of OF text,
%
%   data = y	  \phi, reshaped to vector form.
%   derivFunc	  Function to calculate phi_i^+-.
%   initial	  \phi_0
%
%   delta = ydot  -S(\phi_0)(|\grad \phi| - 1)


% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell 5/27/03
% Calling parameters significantly modified, Ian Mitchell 2/13/04.

  %---------------------------------------------------------------------------
  checkStructureFields(schemeData, 'grid', 'derivFunc', 'initial');

  %---------------------------------------------------------------------------
  grid = schemeData.grid;
  data = reshape(y, grid.shape);

  %---------------------------------------------------------------------------
  % Sign function (smeared) identifies on which side of surface each node lies.
  S = smearedSign(grid, schemeData.initial);

  %---------------------------------------------------------------------------
  % Compute Godunov derivative approximation for each dimension.
  deriv = cell(grid.dim, 1);
  for i = 1 : grid.dim
    [ derivL, derivR ] = feval(schemeData.derivFunc, grid, data, i);

    % For Gudunov's method, check characteristic directions
    %   according to left and right derivative approximations.

    % Both directions agree that flow is to the left.
    flowL = ((S .* derivR <= 0) & (S .* derivL <= 0));

    % Both directions agree that flow is to the right.
    flowR = ((S .* derivR >= 0) & (S .* derivL >= 0));

    % Diverging flow; entropy condition requires choosing deriv = 0
    %   (so we don't actually have to calculate this term).
    %flow0 = ((S .* derivR >  0) & (S .* derivL <  0));

    % Converging flow, need to check which direction arrives first.
    flows = ((S .* derivR <  0) & (S .* derivL >  0));
    if(any(flows(:)))
      conv = find(flows);
      s = zeros(size(flows));
      s(conv) = S(conv) .* (abs(derivR(conv)) - abs(derivL(conv))) ...
                ./ (derivR(conv) - derivL(conv));

      % If s == 0, both directions arrive at the same time.
      %   Assuming continuity, both will produce same result, so pick one.
      flowL(conv) = flowL(conv) | (s(conv) < 0);
      flowR(conv) = flowR(conv) | (s(conv) >= 0);
    end
    
    deriv{i} = derivL .* flowR + derivR .* flowL;
  end


  %---------------------------------------------------------------------------
  % Compute magnitude of gradient.
  mag = zeros(size(grid.xs{1}));
  for i = 1 : grid.dim;
    mag = mag + deriv{i}.^2;
  end
  mag = max(sqrt(mag), eps);

  %---------------------------------------------------------------------------
  % Start with constant term in the reinitialization equation.
  delta = -S;

  % Compute change in function and bound on step size.
  stepBoundInv = 0;
  for i = 1 : grid.dim

    % Effective velocity field (for timestep bounding).
    v = S .* deriv{i} ./ mag;

    % Update just like a velocity field.
    delta = delta + v .* deriv{i};

    % CFL condition using effective velocity.
    stepBoundInv = stepBoundInv + max(abs(v(:))) / grid.dx(i);

  end
  
  %---------------------------------------------------------------------------
  stepBound = 1 / stepBoundInv;
  
  % Reshape output into vector format and negate for RHS of ODE.
  ydot = -delta(:);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function s = smearedSign(grid, data)
% s = smearedSign(grid, data)
%
% Helper function to generated a smeared signum function.
%
% This version is (7.5) in O&F chapter 7.4.

dx = max(grid.dx);
s = data ./ sqrt(data.^2 + dx.^2);
