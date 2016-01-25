% demo1_GLM_temporalStim.m
%
% Demo script for simulating and fitting a single-neuron GLM with  1D
% (temporal) filter with exponential nonlinearity 

% Make sure paths are set (assumes this script called from 'demos' directory)
cd ..; setpaths; cd demos/

%% 1.  Set parameters and display for GLM % =============================

dtStim = .01; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = .001;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
nkt = 30;    % Number of time bins in stimulus filter k
ggsim = makeSimStruct_GLM(nkt,dtStim,dtSp); % Create GLM structure with default params

% === Plot true model params =======================
clf;
ttk = dtStim*(-nkt+1:0)';  % time relative to spike of stim filter taps
subplot(221);plot(ttk, ggsim.k);
title('stimulus kernel'); xlabel('time (s)');

subplot(222); % --------
plot(ggsim.iht, ggsim.iht*0, 'k--', ggsim.iht, ggsim.ih);
title('post-spike kernel h'); axis tight;
xlabel('time after spike (s)');
set(gca, 'ylim',[min(ggsim.ih)*1.1 max(ggsim.ih)*1.5]);

subplot(224); % --------
[iht,ihbasOrthog,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,dtSp);
plot(ggsim.iht, ihbasis);  title('basis for h'); axis tight;


%% 2. Generate some training data  %========================================
slen = 50000; % Stimulus length (frames); more samples gives better fit
swid = 1;  % Stimulus width  (pixels); must match # pixels in stim filter

% Make stimulus
Stim = rand(slen,swid)*2-1;  % Stimulate model to long, unif-random stimulus

% Simulate model
[tsp,sps,Itot,Istm] = simGLM(ggsim,Stim);  % run model

% --- Make plot of first 0.5 seconds of data --------
tlen = 0.5;
ttstim = dtStim:dtStim:tlen; iistim = 1:length(ttstim);
subplot(311); 
plot(ttstim,Stim(iistim));
title('stimulus');

subplot(312);
ttspk = dtSp:dtSp:tlen; iispk = 1:length(ttspk);
spinds = sps(iispk)>0;
semilogy(ttspk,exp(Itot(iispk)),ttspk(spinds), exp(Itot(spinds)), 'ko');
ylabel('spike rate (sp/s)');
title('conditional intensity (and spikes)');

subplot(313); 
Isp = Itot-Istm; % total spike-history filter output
plot(ttspk,Istm(iispk), ttspk,Isp(iispk)); axis tight;
legend('k output', 'h output'); xlabel('time (s)');
ylabel('log intensity'); title('filter outputs');

%% 3. Setup fitting params %===================================================

% Compute the STA
sps_coarse = sum(reshape(sps,[],slen),1)'; % bin spikes in bins the size of stimulus
sta = simpleSTC(Stim,sps_coarse,nkt); % Compute STA
sta = reshape(sta,nkt,[]); % reshape it to match dimensions of true filter

% Set mask (if desired)
exptmask= []; %[1 slen*dtStim];  % data range to use for fitting (in s).

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

%% 4. Do ML fitting %=====================

opts = {'display', 'iter', 'maxiter', 100}; % options for fminunc
[gg1, negloglival] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)

%% 5. Plot results ======================================================

ttk = -nkt+1:0; % time bins for stimulus filter

subplot(221);  % True filter 
plot(ttk, ggsim.k, 'k', ttk, sta./norm(sta)*norm(ggsim.k), ttk, gg1.k, 'r');
title('Stim filters');
legend('k_{true}', 'k_{STA}', 'k_{ML}', 'location', 'northwest');
% ----------------------------------
subplot(222); 
plot(ggsim.iht, ggsim.ih, gg1.iht, gg1.ih,'r', ggsim.iht, ggsim.iht*0, 'k--');
title('post-spike kernel');
axis tight;
legend('h_{true}', 'h_{ML}', 'location', 'southeast');
% ----------------------------------
subplot(224); 
plot(ggsim.iht, exp(ggsim.ih), gg1.iht, exp(gg1.ih),'r', ggsim.iht,ggsim.iht*0+1,'k--');
title('exponentiated post-spike kernels');
xlabel('time since spike (s)');
ylabel('gain'); axis tight;
legend('h_{true}', 'h_{ML}', 'location', 'southeast');

% Errors in STA and ML estimate (subspace angle between true k and estimate)
fprintf('Filter estimation error (in radians)\n  sta: %.3f\n   ML: %.3f\n', ...
    subspace(ggsim.k,sta), subspace(ggsim.k,gg1.k));
