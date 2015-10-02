function DerivCheck(funptr, X0, opts, varargin);
% DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
%
%  Checks the derivative of a function 'funptr' at a point X0, 
%  for purposes of optimization
%
%  Call with same arguments as you would call the optimization routine.

args = {funptr, X0, varargin{:}};  % Arguments to function
[f, JJ] = feval(args{:});  % Evaluate function at X0

tol = 1e-6;  % Size of numerical step to take
rr = randn(length(X0),1)*tol;  % Generate small random-direction vector

X2 = X0+rr;  % Second point at which to evaluate function
args2 = {funptr, X2, varargin{:}};
[f2, JJ2] = feval(args2{:});  % Evaluate at nearby point

fprintf('Derivs: Analytic vs. Finite Diff = [%.4e, %.4e]\n', dot(rr, JJ), f2-f);
