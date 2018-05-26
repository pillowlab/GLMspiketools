function [fun,pp,neglogli,splPrs,Mspline] = fit_tspline_poisson(x,y,ss,dtbin,prs0)
% [fun,pp,neglogli,splPrs,Mspline] = fit_tspline_poisson(x,y,ss,dtbin,prs0)
% 
%  Fit a function y = f(x) with a cubic spline, defined using a set of
%  breaks, smoothness and extrapolation criteria, by maximizing Poisson
%  likelihood:  y ~ Poiss(g(spline(x))).
%
% Inputs:
%       x [Nx1] - input variable
%       y [Nx1] - output (count) variable
%   ss [struct] - spline struct with fields:
%        .breaks - breaks between polynomial pieces
%        .smoothness - derivs of 1 less are continuous
%           (e.g., smoothness=3 -> 2nd derivs are continuous)
%        .extrapDeg - degree polynomial for extrapolation on ends
%        .tfun - nonlinear transfer function (forcing positive outputs)
%   dtbin [1x1]- time bin size (OPTIONAL; assumed 1 otherwise)
%   prs0 [Mx1] - initial guess at spline params (OPTIONAL)
% 
% Outputs:
%   fun - function handle for nonlinearity (uses 'ppval')
%   pp - piecewise polynomial structure
%   neglogli - negative loglikelihood at parameter estimate
%   splinePrs -  vector x such that spline coeffs = Mspline*x
%   Mspline - matrix for converting params to spline coeffs
%
% Note: the fitted nonlinearity can be evaluated via:
%       y = ss.tfun(ppval(pp,xx));
%
% last updated: 26/03/2012 (JW Pillow)

if nargin<4
    dtbin=1;
end

sortflag=1;  % allow sorted design matrix (slightly faster)

% Compute design matrix
[Xdesign,Ysrt,Mspline] = mksplineDesignMat(x,y,ss,sortflag);
nprs = size(Xdesign,2);

if nargin<5
    prs0 = randn(nprs,1)*.1;
end

% Set up loss function
floss = @(prs)neglogli_tspline_poiss(prs,Xdesign,Ysrt,ss.tfun,dtbin);
% HessCheck(floss,prs0); % Check that analytic grad and Hessian are correct

% minimize negative log-likelihood using Newton's method
[splPrs,neglogli] = fminNewton(floss,prs0);  

% Create function handle for resulting tspline
[~,pp] = mksplinefun(ss.breaks, Mspline*splPrs); 
splfun = @(x)ppfun(pp,x);
tfun = ss.tfun;
fun = @(x)tsplinefun(x,splfun,tfun);


% % -----------------------------------------
% % If desired, use Matlab's fminunc instead:
% opts = optimset('display','off','gradobj','on','Hessian','off','maxiter',5000,'maxfunevals',5000,'largescale','off');
% splinePrs2 = fminunc(floss,prs0,opts);
% [floss(splinePrs) floss(splinePrs2)]  % Compare results

