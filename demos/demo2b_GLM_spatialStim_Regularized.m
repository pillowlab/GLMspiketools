% demo2b_GLM_spatialStim_Regularized.m
%
% Test code for simulating and fitting the GLM with a 2D stimulus filter
% (time x 1D space), regularized with a Gaussian prior.

% Make sure paths are set (assumes this script called from 'demos' directory)
cd ..; setpaths_GLMspiketools; cd demos/


%% 1.  Set parameters and display for GLM ============= % 

dtStim = .01; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = .001;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);

% Make a temporal filter
nkt = 20;    % Number of time bins in stimulus filter k
kt = (normpdf(1:nkt,3*nkt/4,1.5)-.5*normpdf(1:nkt,nkt/2,3))';

% Make a spatial filter;
nkx = 10; 
xxk = (1:1:nkx)'; % pixel locations
ttk = dtStim*(-nkt+1:0)'; % time bins for filter
kx = 1./sqrt(2*pi*4).*exp(-(xxk-nkx/2).^2/5);
ktrue = kt*kx'; % Make space-time separable filter
ktrue = ktrue./norm(ktrue(:));

% Insert into glm structure (created with default history filter)
ggsim = makeSimStruct_GLM(nkt,dtStim,dtSp); % Create GLM structure with default params
ggsim.k = ktrue; % Insert into simulation struct
ggsim.dc = 3; 

% === Make Fig: model params =======================
clf; subplot(3,3,[1 4]); % ------------------------------------------
plot(kt,ttk);  axis tight;
set(gca, 'ydir', 'reverse'); ylabel('time (frames)');

subplot(3,3,[2 3 5 6]); % --------------------------------------
imagesc(xxk,ttk,ggsim.k); 
colormap gray; axis off; 
title('stimulus kernel k');

subplot(3,3,8:9); % ----------------------------------------
plot(xxk,kx); axis tight;
set(gca, 'xlim', [.5 nkx+.5]); xlabel('space (pixels)');


%% 2. Generate some training data ========================================

% generate stimulus
slen = 10000; % Stimulus length (frames);  More samples gives better fit
gfilt = normpdf(-3:3, 0, .8); 
Stim = conv2(randn(slen,nkx),gfilt,'same'); %  convolve in space
Stim = conv2(Stim,gfilt','same'); % convolve in time

[tsp,sps,Itot,Istm] = simGLM(ggsim,Stim);  % run model
nsp = length(tsp);
rlen = length(sps);

% --- Make plot of first 0.5 seconds of data --------
tlen = 0.5;
ttstim = dtStim:dtStim:tlen; iistim = 1:length(ttstim);
ttspk = dtSp:dtSp:tlen; iispk = 1:length(ttspk);
spinds = sps(iispk)>0;
subplot(311); 
imagesc([0 tlen], [1 nkx], Stim(iistim,:)'); 
title('stimulus'); ylabel('pixel');

subplot(312);
semilogy(ttspk,exp(Itot(iispk)),ttspk(spinds), exp(Itot(spinds)), 'ko');
ylabel('spike rate (sp/s)');
title('conditional intensity (and spikes)');

subplot(313); 
Isp = Itot-Istm; % total spike-history filter output
plot(ttspk,Istm(iispk), ttspk,Isp(iispk)); axis tight;
legend('k output', 'h output'); xlabel('time (s)');
ylabel('log intensity'); title('filter outputs');



%% 3. Fit GLM with max likelihood

% Compute the STA
sps_coarse = sum(reshape(sps,[],slen),1)'; % bin spikes in bins the size of stimulus
sta0 = simpleSTC(Stim,sps_coarse,nkt); % Compute STA
sta = reshape(sta0,nkt,[]); % reshape it to match dimensions of true filter

exptmask= [];  % Not currently supported!

%  Initialize params for fitting --------------
% Set params for fitting, including bases 
nkbasis = 8;  % number of basis vectors for representing k
nhbasis = 8;  % number of basis vectors for representing h
hpeakFinal = .1;   % time of peak of last basis vector for h
gg0 = makeFittingStruct_GLM(dtStim,dtSp,nkt,nkbasis,sta,nhbasis,hpeakFinal);
gg0.sps = sps;  % Insert binned spike train into fitting struct
gg0.mask = exptmask; % insert mask (optional)
gg0.ihw = randn(size(gg0.ihw))*1; % initialize spike-history weights randomly

% Compute conditional intensity at initial parameters 
[negloglival0,rr] = neglogli_GLM(gg0,Stim);
fprintf('Initial negative log-likelihood: %.5f\n', negloglival0);

% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 1000};
[gg1, negloglival1] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% 4. Fit GLM with iid Gaussian ("ridge regression") prior.

% Divide data into training and test
trainfrac = 0.8;
ntrain = ceil(trainfrac*slen);
ntrainsps = ceil(trainfrac*rlen);
stimtrain = Stim(1:ntrain,:); stimtest = Stim(ntrain+1:end,:);
spstrain = sps(1:ntrainsps); spstest = sps(ntrainsps+1:end); 

% Set prior covariance matrix ---------------
nparams = numel(gg0.kt);
Ceye = speye(nparams);

% Make grid of lambda (ridge parameter) values -----
lamvals = 2.^(-1:10);
nlam = length(lamvals);
testlogli = zeros(nlam,1);
trainlogli = zeros(nlam,1);
ggridge = cell(nlam,1);

% Initialize params --------
gg0_train = gg0;
gg0_train.sps = spstrain;

% Find MAP estimate for each value of ridge parameter
fprintf('====Doing grid search for ridge parameter====\n');
for jj = 1:nlam
    fprintf('lambda = %.1f\n',lamvals(jj));
    Cinv = lamvals(jj)*Ceye;  % inverse covariance of prior
    
    % Do MAP estimation of model params
    [ggridge{jj}, negloglival] = MAPfit_GLM(gg0_train,stimtrain,Cinv,opts); % do ML (requires optimization toolbox)
    trainlogli(jj) = -negloglival;
    
    % Compute test log-likelihood
    gg2tst = ggridge{jj};
    gg2tst.sps = spstest;
    testlogli(jj) = -neglogli_GLM(gg2tst,stimtest);
end

% Plot test log likelihood
clf; semilogx(lamvals,testlogli, '-o');
xlabel('lambda'); ylabel('test log-likelihood'); title('ridge prior');

% Select best ridge param and fit on all data.
[~,imax] = max(testlogli);
[gg2,neglogli2] = MAPfit_GLM(gg1,Stim,lamvals(imax)*Cinv,opts);

%% 5. Fit GLM with smoothing prior.

% Make matrices to compute column and row pairwise differences
Dx = spdiags(ones(nkx,1)*[-1 1],0:1,nkx-1,nkx);
Dt = spdiags(ones(nkt,1)*[-1 1],0:1,nkt-1,nkt);
Dtbas = Dt*gg1.ktbas;
Mt = kron(speye(nkx),Dtbas'*Dtbas);  % computes squared diffs along columns
Mx = kron(Dx,gg1.ktbas); % computes squared diffs along row
Mx = Mx'*Mx;
C0 = Mt + Mx;  % default covariance matrix

% Make grid of lambda (smoothing parameter) values -----
lamvals_sm = 2.^(0:13);
nlam = length(lamvals_sm);
testlogli_sm = zeros(nlam,1);
trainlogli_sm = zeros(nlam,1);
ggsmooth = cell(nlam,1);

% Find MAP estimate for each value of ridge parameter
for jj = 1:nlam
    Cinv = lamvals_sm(jj)*C0; % inverse covariance of prior
    
    % Do MAP estimation of model params
    [ggsmooth{jj}, negloglival] = MAPfit_GLM(gg0_train,stimtrain,Cinv,opts); % do ML (requires optimization toolbox)
    trainlogli_sm(jj) = -negloglival;
    
    % Compute test log-likelihood
    gg3tst = ggsmooth{jj};
    gg3tst.sps = spstest;
    testlogli_sm(jj) = -neglogli_GLM(gg3tst,stimtest);
end

% Plot test log likelihood
clf; semilogx(lamvals_sm,testlogli_sm, '-o');
xlabel('lambda'); ylabel('test log-likelihood'); title('smoothing prior');

% Select best ridge param and fit on all data.
[~,imax] = max(testlogli_sm);
[gg3,neglogli3] = MAPfit_GLM(gg2,Stim,lamvals_sm(imax)*C0,opts);


%% 6. Plot results ====================

subplot(231);  % True filter  % ---------------
imagesc(ggsim.k); colormap gray;
title('True Filter');ylabel('time');
subplot(232);  % sta % ------------------------
imagesc(gg1.k);
title('ML estimate'); ylabel('time');

subplot(234); % sta-projection % ---------------
imagesc(gg2.k); title('MAP: ridge prior'); xlabel('space'); ylabel('time');

subplot(235); % estimated filter % ---------------
imagesc(gg3.k); title('MAP: smoothing prior'); xlabel('space'); 

subplot(236); % ----------------------------------
plot(ggsim.iht,exp(ggsim.ih),'k--', gg1.iht,exp(gg1.ihbas*gg1.ihw),...
    gg2.iht, exp(gg2.ihbas*gg2.ihw), ...
    gg3.iht, exp(gg3.ihbas*gg3.ihw) ...
    ); axis tight;
title('post-spike kernel');  xlabel('time after spike (s)');
legend('true','ML','ridge','smooth');

r2fun = @(k)(1-sum((k(:)-ktrue(:)).^2)./sum(ktrue(:).^2)); % error function
fprintf('\nK filter R^2:\n ----------\n ML = %9.3f\n Ridge = %6.3f\n Smooth = %4.3f\n',...
    r2fun(gg1.k), r2fun(gg2.k),r2fun(gg3.k));
