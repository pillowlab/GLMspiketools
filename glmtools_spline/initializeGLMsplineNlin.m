function [ggnew,negloglival] = initializeGLMsplineNlin(gg,Stim,prs)
% [ggnew,negloglival] = initializeGLMsplineNlin(gg,Stim,prs)
% 
% Initializes parameters for a spline nonlinearity using GLM log-likelihood 
%
% Inputs: 
%    gg = param struct
%    Stim = stimulus
%    spline_prs = structure parameters for spline nonlinearity
%        ---- Fields: ----
%        nknots = # of knots,
%        epprob = [low,high] quantiles for first and last knots
%        smoothness = 1+degree of smoothness (3 -> 2nd order continuous)
%        extrapDeg = [lo,hi] degree of polynomials for extrapolation outside breaks
%      
%  Outputs:
%     ggnew = new param struct (with ML params);
%     negloglival = negative log-likelihood at ML estimate

% Compute net current output from GLM
[~,~,~,Itot] = neglogli_GLM(gg,Stim);
Itot = Itot-gg.dc; % remove DC from injected current

% Set knots and insert piece-wise linear nonlinearity.
knotrange  = quantile(Itot,prs.epprob);
knots = linspace(knotrange(1),knotrange(2), prs.nknots);
prebreaks = ((prs.extrapDeg(1):1)-2)/10; % extra left breaks (spaced 0.1 apart)
postbreaks = (1:(2-prs.extrapDeg(2)))/10; % extra right breaks (spaced 0.1 apart)
breaks = [prebreaks+knots(1), knots, knots(end)+postbreaks];

% make spline structure
ss = struct('breaks',breaks, 'smoothness',prs.smoothness, 'extrapDeg',prs.extrapDeg,'tfun',@logexp1);
[fspline,pp,negloglival] = fit_tspline_poisson(Itot,gg.sps,ss,gg.dtSp);

% Insert params into glm struct
ggnew = gg; % initialize
ggnew.dc = 0; % set DC term to zero
ggnew.ppstruct = pp;  % piecewise polynomial parameters
ggnew.splineprs = ss; % spline hyperparameters
ggnew.nlfun = fspline; % insert nonlinearity

% % Optional: make plot comparing old and new nonlinearity
% xx = linspace(min(II),max(II),100);
% plot(xx,fspline(xx),breaks,fspline(breaks),'o',xx,gg.nlfun(xx));
