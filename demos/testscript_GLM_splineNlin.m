% testscript_GLM_splineNlin.m
%
% Test code for simulating and fitting a GLM with a nonlinearity
% parametrized by a cubic spline (instead of a fixed nonlinearity)
%
% Code blocks:
%  1. Create GLM with non-exponential nonlinearity
%  2. Fit regular GLM
%  3. Fit GLM with spline-parametrized nonlinearity
%  4. Generate repeated responses to a short stimulus, plot rasters


global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100; 

%% 1.  Set up GLM with arbitrary nonlinearity & make training data
DTsim = .01; % Bin size for simulating model & computing likelihood.
nkt = 15;  % Number of time bins in filter;
nkx = 5;   % Number of spatial pixels
ggsim = makeSimStruct_GLM(nkt,DTsim); % Create GLM structure with default params
ggsim.dc = 0;
kt = ggsim.k;

% Create spatial filter
kx = (normpdf(1:nkx, (nkx+1)/2,1)*2 - normpdf(1:nkx, (nkx+1)/2,2))';
k = 5*kt*kx';  % Full filter is space-time separable

% Set up post-spike kernel
k_rank = 1;
ggsim = makeFittingStruct_GLMbi(k,k_rank,DTsim);
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,DTsim);
ih = ihbas'*(ihbasis*[-2;-2;0;.55;-.25]);
%plot(iht, ihbas*ih);
ggsim.ih = ih;

% Create nonlinearity
ggsim.nlfun = @(x)((100*abs(x/2+.3).^4./(1+abs(x/2+.3).^4)+2).*(1-.8.*(x<0)));

% Simulate response to a long stimulus
slen = 1000; % Stimulus length (frames);  More samples gives better fit
Stim = randn(slen,nkx);  %  Run model on long, binary stimulus
[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model
nsp = length(tsp)
ggsim.tsp = tsp; 
ggsim.tspi = 1;

% Compute STA 
sta0 = simpleSTC(Stim,tsp,nkt);
sta = reshape(sta0,nkt,[]);

% make plots --------
figure(1);  
subplot(221); imagesc(ggsim.k); title('True filter');
[lv0,rr0,tt,Itot0] = neglogli_GLM(ggsim,Stim);
xnlin0 = min(Itot0):.01:max(Itot0);
subplot(223); imagesc(sta); title('STA');
subplot(222); plot(xnlin0, ggsim.nlfun(xnlin0)); title('nonlinearity');
axis([-12 12 0 180])
subplot(224); plot(ggsim.iht, ggsim.ihbas*ggsim.ih); title('post-spike h');
axis tight;
colormap gray;
drawnow;


%% 2.  Fit GLM with fixed nonlinearity (eg, exponential or log(1+exp(x)) ) 
% Initialize by fitting exponential GLM;
gg0 = makeFittingStruct_GLMbi(sta,k_rank,DTsim);
%gg0.nlfun = @expfun;  % Exponential nonlinearity
gg0.nlfun = @logexp3  % f(x) = log(1+exp(x)).^3
gg0.tsp = tsp; 
gg0.tspi = 1;

% ---- Compare initial and true log-likelihood, if desired ----------
% [logliTrue, rrT] = neglogli_GLM(ggsim,Stim);
% [logli0, rr0,tt0] = neglogli_GLM(gg0,Stim);
% fprintf('Logli: trueparams=%.3f, Initial=%.3f\n', logliTrue,logli0);
% % -----------

% Do ML estimation of model params with standard (bilinear-filter) GLM
opts = {'display', 'iter', 'maxiter', 100};
[gg1, loglival] = MLfit_GLMbi(gg0,Stim,opts); 

% ------ make plots ---------------------------
figure(2);  clf;
subplot(221); 
imagesc(ggsim.k); title('true filter');
subplot(223); 
imagesc(gg1.k); title('filter estimate (regular GLM)');

% Compute x values relevant for nonlinearity with a call to "neglogli"
[lv1,rr1,tt,Itot1] = neglogli_GLM(gg1,Stim);
xnlin1 = min(Itot1):.01:max(Itot1);

subplot(222); 
plot(xnlin1, gg1.nlfun(xnlin1), 'b'); title('nonlinearity');

subplot(224); 
plot(ggsim.iht, normalizecols(ggsim.ihbas*ggsim.ih), 'k--', ...
    gg1.iht, normalizecols(gg1.ihbas*gg1.ih));
axis tight; 
title('post-spike h (normalized)');
drawnow;
colormap gray;

%% 3.  Fit GLM with spline-parametrized nonlinearity

% Spline parameters
nknots = 5;  % number of knots (# cubic segments - 1)
epprob = [1e-3, 1-1e-3];  % knot endpoints (in terms of cum distribution of filter output)
smthness = 3;  % polynomial deg+1 for spline continuity (3 => 2nd order continuity)
extrapDeg=[1 1]; % degree of polynomials on each end (for extrapolation)
minval = 1e-10; % Min value of spline nonlinearity (rectification)

splinePrs.nknots = nknots;  % total # knots is this + 1
splinePrs.epprob = epprob;
splinePrs.smthness = smthness;
splinePrs.extrapDeg = extrapDeg;
splinePrs.minval = minval;

% Initialize spline nonlinearity (fit with filter params fixed)
opts = {'maxiter', 20};
gg2 = initializeSplineNlin_GLM(gg1,Stim,splinePrs,opts);

% ----------------------
% Iteratively fit filter and spline params
npasses = 10;  % # of passes of iterative fitting
opts = {'maxiter', 20, 'display', 'iter'};
gg2 = MLfit_GLMspline(gg2,Stim,opts,npasses);
% Run this part again if dLoss isn't < 0.1

% ---------- make plots ----------------------------------------------
figure(3);  clf;
subplot(221); 
imagesc(ggsim.k); title('true filter');
subplot(223); 
imagesc(gg2.k); title('filter estimate (GLM w/ spline)');

% Compute x values relevant for nonlinearity with a call to "neglogli"
[lv2,rr2,tt,Itot2] = neglogli_GLM(gg2,Stim);
xnlin2 = min(Itot2):.01:max(Itot2);
knots = gg2.ppstruct.breaks;

subplot(222); 
plot(xnlin2, gg2.nlfun(xnlin2),knots,gg2.nlfun(knots), 'ro'); 
title('nonlinearity');

subplot(224); 
plot(ggsim.iht, normalizecols((ggsim.ihbas*ggsim.ih)), 'k--', ...
    gg2.iht, normalizecols((gg2.ihbas*gg2.ih)));
axis tight; 
title('post-spike h (normalized)');
colormap gray;

Estim_Error = [subspace(ggsim.k(:),sta(:)), ...
    subspace(ggsim.k(:),gg1.k(:)), ...
    subspace(ggsim.k(:),gg2.k(:))]
% Note that the magnitude of the estimated filters may differ from the 
% true filters, because nonlinearity can interacts with filter amplitude.


%% 4.  Simulate responses to a short repeated stimulus
slen = 100;  % Stimulus length
nrpts = 50;  % # repeats
Stim = randn(slen,nkx);
tspsim = []; tsp1 = []; tsp2 = [];

fprintf('Simulating: Rpt #\n');
for j = 1:nrpts;
    fprintf(' %d', j);
    tspsim{j} = simGLM(ggsim,Stim);
    tsp1{j} = simGLM(gg1,Stim);
    tsp2{j} = simGLM(gg2,Stim);
end
fprintf('\n');

%% --- Make figure ---------------------
figure(4); 
clf;
plotraster(tspsim,tsp1,tsp2,[0 slen]);
hold on;
plot([0 slen], nrpts*[1 1], 'k', ...
    [0 slen], nrpts*[2 2], 'k');
hold off;
set(gca, 'ytick', nrpts/2 + nrpts*[0:2], ...
    'yticklabel', {'true', 'GLM', 'spline GLM'});
xlabel('time (frames)')