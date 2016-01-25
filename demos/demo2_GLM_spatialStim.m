% demo2_GLM_spatialStim.m
%
% Test code for simulating and fitting the GLM with a 2D stimulus filter
% (time x 1D space), with both traditional and bilinear parametrization
% of the stimulus kernel.

% Make sure paths are set (assumes this script called from 'demos' directory)
cd ..; setpaths; cd demos/


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
k = kt*kx'; % Make space-time separable filter
k = k./norm(k(:))/1.5;

% Insert into glm structure (created with default history filter)
ggsim = makeSimStruct_GLM(nkt,dtStim,dtSp); % Create GLM structure with default params
ggsim.k = k; % Insert into simulation struct
ggsim.dc = 2; 

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

slen = 10000; % Stimulus length (frames);  More samples gives better fit
Stim = round(rand(slen,nkx))*2-1;  %  Run model on long, binary stimulus
[tsp,sps,Itot,Istm] = simGLM(ggsim,Stim);  % run model
nsp = length(tsp);

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



%% 3. Fit GLM (traditional version) via max likelihood

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
opts = {'display', 'iter', 'maxiter', 100};
[gg1, negloglival1a] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% 5. Fit GLM ("bilinear stim filter version") via max likelihood

%  Initialize params for fitting --------------
k_rank = 1; % Number of column/row vector pairs to use
gg0b = makeFittingStruct_GLMbi(k_rank,dtStim,dtSp,nkt,nkbasis,sta,nhbasis,hpeakFinal);
gg0b.sps = sps;
gg0b.mask = exptmask;
logli0b = neglogli_GLM(gg0b,Stim); % Compute logli of initial params
fprintf('Initial value of negative log-li (GLMbi): %.3f\n', logli0b);

% Do ML estimation of model params
opts = {'display', 'iter'};
[gg2, negloglival2] = MLfit_GLMbi(gg0b,Stim,opts); % do ML (requires optimization toolbox)


%% 6. Plot results ====================
%figure(3);

subplot(231);  % True filter  % ---------------
imagesc(ggsim.k); colormap gray;
title('True Filter');ylabel('time');
subplot(232);  % sta % ------------------------
imagesc(sta);
title('raw STA'); ylabel('time');

subplot(233); % sta-projection % ---------------
imagesc(gg0.k); title('low-rank STA');

subplot(234); % estimated filter % ---------------
imagesc(gg1.k); title('ML estimate: full filter'); xlabel('space'); ylabel('time');

subplot(235); % estimated filter % ---------------
imagesc(gg2.k); title('ML estimate: bilinear filter'); xlabel('space'); 

subplot(236); % ----------------------------------
plot(ggsim.iht,exp(ggsim.ih),'k', gg1.iht,exp(gg1.ihbas*gg1.ihw),'b',...
    gg2.iht, exp(gg2.ihbas*gg2.ihw), 'r'); axis tight;
title('post-spike kernel');  xlabel('time after spike (s)');
legend('true','GLM','bilinear GLM');

% Errors in STA and ML estimate
ktmu = normalizecols([mean(ggsim.k,2),mean(gg1.k,2),mean(gg2.k,2)]);
kxmu = normalizecols([mean(ggsim.k)',mean(gg1.k)',mean(gg2.k)']);

msefun = @(k)(sum((k(:)-ggsim.k(:)).^2)); % error function
fprintf(['K-filter errors (GLM vs. GLMbilinear):\n', ...
    'Temporal error:  %.3f    %.3f\n', ...
    ' Spatial error:  %.3f    %.3f\n', ...
    '   Total error:  %.3f    %.3f\n'], ...    
    subspace(ktmu(:,1),ktmu(:,2)), subspace(ktmu(:,1),ktmu(:,3)), ...
    subspace(kxmu(:,1),kxmu(:,2)), subspace(kxmu(:,1),kxmu(:,3)), ...
    msefun(gg1.k), msefun(gg2.k));
