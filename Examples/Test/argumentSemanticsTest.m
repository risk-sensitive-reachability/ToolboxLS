function argumentSemanticsTest(loops, matSize)
% argumentSemanticsTest: test Matlab's argument passing speed.
%
%   argumentSemanticsTest(loops, matSize)
%
% Script file to test the effectiveness of Matlab's purported pass by value
%   semantics with pass by reference speed.
%
% Specifically, Matlab uses pass by value semantics, but the documentation
%   claims that no copy of arguments are made (ie pass by reference) until
%   such time as the argument is modified, at which point a copy is made.
%
% A related question is whether statements like B = A and A = A(:)
%   cause copies to be made immediately.
%
% My testing seems to indicate that Matlab lives up to its billing:
%   the functions copyAndModify and Modify below both take much 
%   longer than the functions noOp, makeVectorCopy, 
%   makeVector and reshape2Column.
%
% parameters:
%   loops        Number of loops to execute in order to get timings.
%   matSize      Size of matrices to use in the timing 
%                  (should be an even integer).

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/6/04

%---------------------------------------------------------------------------
  % File handle for the screen.
  screen = 1;
  fprintf(screen, 'Making repeating each function call %d times\n', loops);
  
  functions = { @noOp; @copyAndModify; @modify; @makeVectorCopy; ...
                @makeVector; @reshape2Column };
  
  A = rand(matSize);
  
  for i = 1 : length(functions)
    startTime = cputime;
    for j = 1 : loops
      B = feval(functions{i}, A);
    end
    endTime = cputime;
    fprintf(screen, '  Function %s finished in %g seconds\n', ...
            func2str(functions{i}), endTime - startTime);
  end
  
  
  
  
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function B = noOp(A)
% Basic function that performs no useful operations.
  
  B = A;
  
    
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function B = copyAndModify(A)
% Copies the input argument, and modifies the copy.
  
  B = A;
  B(1,1) = 1;
    
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function B = modify(A)
% Modifies the input argument, but makes no copy.
  
  A(1,1) = 2;
  B = A;
  
    
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function B = makeVectorCopy(A)
% Copy is argument reshaped into a vector with the (:) notation.
  
  B = A(:);
  
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function B = makeVector(A)
% Reshapes the argument into a vector with the (:) notation.
  
  A = A(:);
  B = A(1);
  
    
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function B = reshape2Column(A)
% Reshapes the argument into two columns.
  
  lengthA = prod(size(A));
  A = reshape(A, lengthA / 2, 2);
  B = A(1,1);
  
  
