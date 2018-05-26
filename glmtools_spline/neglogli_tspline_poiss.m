function [L,dL,ddL] = neglogli_tspline_poiss(prs,X,Y,g,dtbin)
% [L,dL,ddL] = neglogli_tspline_poiss(prs,X,Y,g,dtbin)
%
% Compute negative log-likelihood of data Y given rate g(X*prs).
%
% INPUT:
%       prs [M x 1] - parameters 
%         X [N x M] - design matrix
%         Y [N x 1] - observed Poisson random variables
%                 g - handle for transfer function
%     dtbin [1 x 1] - size of time bins of Y 
%
% OUTPUT: 
%         L [1 x 1] - negative log-likelihood
%        dL [M x 1] - gradient
%       ddL [M x M] - Hessian (2nd deriv matrix)
%
% Last modified: 25 May 2018 (JW Pillow)


% Project parameters onto design matrix
z = X*prs;
etol = 1e-100;

if nargout==1
    % Compute neglogli
    f = g(z);
    f(f<etol)=etol;
    L = -Y'*log(f) + sum(f)*dtbin;
elseif nargout == 2
    % Compute neglogli & Gradient
    [f,df] = g(z);
    f(f<etol)=etol;
    L = -Y'*log(f) + sum(f)*dtbin;
    % grad
    wts = (dtbin*df-(Y.*df./f));
    dL = X'*wts;
elseif nargout == 3
    % Compute neglogli, Gradient & Hessian
    [f,df,ddf] = g(z);
    f(f<etol)=etol;
    L = -Y'*log(f) + sum(f)*dtbin;
    % grad
    wts = (dtbin*df-(Y.*df./f));
    dL = X'*wts;
    ww = dtbin*ddf-Y.*(ddf./f-(df./f).^2);
    ddL = X'*bsxfun(@times,X,ww);
end
