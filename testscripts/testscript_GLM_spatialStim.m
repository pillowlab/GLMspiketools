% testscript_GLMbi.m
%
% Test code for simulating and fitting the GLM with a 2D stimulus filter
% (1D time by 1D spatial), with both traditional and bilinear parametrization
% of the stimulus kernel.
%
% Code Blocks:  
%   1. Set up model params and plot
%   2. Show samples from simulated model
%   3. Generate some training data 
%   4. Fit traditional GLM via maximum-likelihood (ML)
%   5. Fit GLM with bilinearly-parametrized filter ("GLMbi") via ML

global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100; 

%% 1.  Set parameters and display for GLM ============= % 

DTsim = .01; % Bin size for simulating model & computing likelihood.
nkt = 15;  % Number of time bins in filter;
ttk = [-nkt+1:0]';
ggsim = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
kt = ggsim.k;  % Temporal filter

% Make a spatial filter;
nkx = 10; 
xxk = [1:1:nkx]';
kx = 1./sqrt(2*pi*4).*exp(-(xxk-nkx/2).^2/5);
Filt = kt*kx'; % Make space-time separable filter
ggsim.k = Filt./norm(Filt(:))*3; % Insert into simulation struct

figure(1);  % === Make Fig: model params =======================
subplot(3,5,[1,6]); % ------------------------------------------
plot(kt,ttk);  axis tight;
set(gca, 'ydir', 'reverse');
ylabel('time (frames)');

subplot(3,5,[2,3,7,8]); % --------------------------------------
imagesc(xxk,ttk,ggsim.k); 
colormap gray;
axis off; 
title('stimulus kernel k');

subplot(3,5,[12,13]); % ----------------------------------------
plot(xxk,kx); axis tight;
set(gca, 'xlim', [.5 nkx+.5]);
xlabel('space (pixels)');

subplot(3,5,4:5); % --------------------------------------------
plot(ggsim.iht, ggsim.iht*0, 'k--', ggsim.iht, ggsim.ih); 
title('post-spike kernel h');
set(gca, 'xlim', ggsim.iht([1 end]),...
    'ylim',[min(ggsim.ih)*1.1 max(ggsim.ih)*1.5]);
subplot(3,5,9:10); % -------------------------------------------
[iht,ihbasOrthog,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,DTsim);
plot(ggsim.iht, ihbasis);
axis tight;
xlabel('time after spike (frames)');
title('basis for h');

%% 2. Make GWN stimulus & simulate the glm model response. ========= %

slen = 50; % Stimulus length (frames) & width (# pixels)
swid = size(ggsim.k,2);
Stim = randn(slen,swid);  % Gaussian white noise stimulus
[tsp, vmem,Ispk] = simGLM(ggsim, Stim);

% ==== Make Figure ========
figure(2); 
tt = [DTsim:DTsim:slen]';
subplot(221); %------------------------
imagesc(Stim'); 
colormap gray; axis image;
title('GWN stimulus');
xlabel('time');
ylabel('space');

subplot(223); %------------------------
plot(tt, vmem-Ispk, tt, Ispk, 'r', tsp, .1*ones(size(tsp)), 'r.');
title('stim- and spike-induced currents'); axis tight;
xlabel('time (frames)');
subplot(222); %------------------------
plot(tt, vmem, tsp, max(vmem)*ones(size(tsp)), 'ro');
title('net voltage (dots = spike times)'); 
axis tight;
% -----Run repeat simulations ------
nrpts = 5;        % number of repeats to draw
subplot(224);
for j = 1:nrpts;
    [tsp1,vmem1] = simGLM(ggsim,Stim);
    plot(tt,vmem1, 'color', rand(3,1)); hold on;
end
axis tight; hold off;
title('repeat responses to same stim');
xlabel('time (frames)');

%% 3. Generate some training data ========================================

slen = 1000; % Stimulus length (frames);  More samples gives better fit
Stim = round(rand(slen,swid))*2-1;  %  Run model on long, binary stimulus
[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model
nsp = length(tsp);

% Compute STA and use as initial guess for k
sta0 = simpleSTC(Stim,tsp,nkt);
sta = reshape(sta0,nkt,[]);

% % ---------------
% % Make param object with "true" params;
% % ---------------
% Filter_rank = 1;
% ggTrue = makeFittingStruct_GLMbi(ggsim.k,Filter_rank,DTsim,ggsim);
% % Set true kernel in this struct (i.e. not represented by default basis).
% [u,s,v] = svd(ggsim.k);  
% ggTrue.k = ggsim.k;
% ggTrue.ktbas = eye(nkt);
% ggTrue.kt = u(:,1);  
% ggTrue.kx = v(:,1)*s(1,1);
% % Insert spike times
% ggTrue.tsp = tsp;
% ggTrue.tspi = 1;  % Start computing likelihood from 1st spike
% 
% % ---------------
% % Check that conditional intensity calc is correct 
% % (if desired, compare to vmem returned by simGLM above)
% [logliTrue, rrTrue,tt] = neglogli_GLM(ggTrue,Stim);
% % ---------------


%% 4. Fit GLM (traditional version) via max likelihood

%  Initialize params for fitting --------------
Filter_rank = 1;
gg0 = makeFittingStruct_GLM(sta,DTsim);
gg0.tsp = tsp;
gg0.tspi = 1;
[logli0,rr0,tt] = neglogli_GLM(gg0,Stim); % Compute logli of initial params (if desired)

% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 100};
[gg1, negloglival] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% 5. Fit GLM ("bilinear stim filter version") via max likelihood

%  Initialize params for fitting --------------
Filter_rank = 1; % Number of column/row vector pairs to use
gg0 = makeFittingStruct_GLMbi(sta,Filter_rank,DTsim);
gg0.tsp = tsp;
gg0.tspi = 1;
[logli0,rr0,tt] = neglogli_GLM(gg0,Stim); % Compute logli of initial params


% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 100};
[gg2, negloglival] = MLfit_GLMbi(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% 6. Plot results ====================
figure(3);

subplot(231);  % True filter  % ---------------
imagesc(ggsim.k); colormap gray;
title('True Filter');ylabel('time');

subplot(232);  % sta % ------------------------
imagesc(sta);
title('raw STA');
ylabel('time');

subplot(233); % sta-projection % ---------------
imagesc(gg0.k)
title('projected STA');

subplot(234); % estimated filter % ---------------
imagesc(gg1.k) 
title('ML estimate: full filter'); xlabel('space'); ylabel('time');

subplot(235); % estimated filter % ---------------
imagesc(gg2.k)
title('ML estimate: bilinear filter'); xlabel('space'); 

subplot(236); % ----------------------------------
plot(ggsim.iht,exp(ggsim.ih),'k', gg1.iht,exp(gg1.ihbas*gg1.ih),'b',...
    gg2.iht, exp(gg2.ihbas*gg2.ih), 'r');
title('post-spike kernel');
axis tight;

% Errors in STA and ML estimate
ktmu = normalizecols([mean(ggsim.k,2),mean(gg1.k,2),mean(gg2.k,2)]);
kxmu = normalizecols([mean(ggsim.k)',mean(gg1.k)',mean(gg2.k)']);
Errs_T = [subspace(ktmu(:,1),ktmu(:,2)), subspace(ktmu(:,1),ktmu(:,3))]
Errs_X = [subspace(kxmu(:,1),kxmu(:,2)), subspace(kxmu(:,1),kxmu(:,3))]

errfun = @(x,y)(sum((x(:)-y(:)).^2));
Errs_Total = [errfun(ggsim.k,gg1.k), errfun(ggsim.k, gg2.k)]
