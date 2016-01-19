% demo1_GLM_temporalStim.m
%
% Test code for simulating and fitting a single-neuron GLM with  1D
% (temporal) filter with exponential nonlinearity 
%
% Code Blocks:  
%   1. Set up model params
%   2. Show simulated model responses to stimuli
%   3. Make training data (for fitting params to data)
%   4. Fit GLM params via maximum-likelihood (requires optimization toolbox)

% Make sure paths are set (assumes this script called from 'demos' directory)
cd ..; setpaths; cd demos/

%% 1.  Set parameters and display for GLM % =============================

dtStim = .01; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = .001;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
nkt = 30;    % Number of time bins in stimulus filter k
ttk = dtStim*(-nkt+1:0)';  % time relative to spike of stim filter taps
ggsim = makeSimStruct_GLM(nkt,dtStim,dtSp); % Create GLM structure with default params

% === Plot true model params =======================
clf;
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


%% 2. Make some training data  %========================================
slen = 5000; % Stimulus length (frames); more samples gives better fit
swid = 1;  % Stimulus width  (pixels); must match # pixels in stim filter

% Make stimulus
Stim = rand(slen,swid)*2-1;  % Stimulate model to long, unif-random stimulus

% Simulate model
[tsp,Itot,Isp] = simGLM(ggsim,Stim);  % run model

% --- Make plot of first 0.5 seconds of data --------
tlen = 0.5;
ttstim = dtStim:dtStim:tlen;
subplot(311); 
plot(ttstim,Stim(1:length(ttstim)), tsp(tsp<tlen), 0, 'ko');
title('stimulus and spike times');

ttspk = dtSp:dtSp:tlen;
iispk = 1:length(ttspk);
subplot(312); 
plot(ttspk,Itot(iispk)-Isp(iispk), ttspk,Isp(iispk)); axis tight;
legend('k output', 'h output');
ylabel('log intensity'); title('filter outputs');

subplot(313);
spinds = round(tsp(tsp<tlen)/dtSp);
semilogy(ttspk,exp(Itot(iispk)),tsp(tsp<tlen), exp(Itot(spinds)), 'ko');
ylabel('spike rate (sp/s)');xlabel('time (s)');
title('conditional intensity');

%% 4. Setup fitting params %===================================================

% Compute STA and use as initial guess for k
sps = binSpTimes(tsp,dtStim,slen*dtStim);  % bin spike times
sta0 = simpleSTC(Stim,sps,nkt); % compute STA
sta = reshape(sta0,nkt,[]);

% Set mask (if desired)
exptmask= [100 slen]*dtStim;  % data range to use for fitting 

% Set params for fitting, including bases 
nkbasis = 8;  % number of basis vectors for representing k
nhbasis = 5;  % number of basis vectors for representing h
hpeak = .1;   % time of peak of last basis vector for h
gg0 = makeFittingStruct_GLM(dtStim,dtSp,nkt,nkbasis,sta*.25,nhbasis,hpeak);
gg0.tsp = tsp;  % Insert spikes into fitting struct
gg0.mask = exptmask;

% Compute conditional intensity at initial parameters 
negloglival0 = neglogli_GLM(gg0,Stim);
fprintf('Initial negative log-likelihood: %.2f\n', negloglival0);

%% 4. Do ML fitting %=====================

opts = {'display', 'iter', 'maxiter', 100}; % options for fminunc
[gg, negloglival] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)

%% 5. Plot results ======================================================

ttk = -nkt+1:0; % time bins for stimulus filter
subplot(221);  % True filter  % ---------------
plot(ttk, ggsim.k, 'k', ttk, sta, ttk, gg.k, 'r');
title('Stim filters');
legend('k_{true}', 'k_{STA}', 'k_{ML}', 'location', 'northwest');

subplot(223);
flts = ([ggsim.k./norm(ggsim.k), sta, gg.k./norm(gg.k)]);
plot(ttk, flts(:,1),'k', ttk,flts(:,2), ttk, flts(:,3), 'r');
title('Normalized Stim filters');
xlabel('time before spike (frames)');

subplot(222); % ----------------------------------
plot(ggsim.iht, ggsim.ih, gg.iht, gg.ih,'r');
title('post-spike kernel');
axis tight;
legend('h_{true}', 'h_{ML}', 'location', 'southeast');


subplot(224); % ----------------------------------
plot(ggsim.iht, exp(ggsim.ih), gg.iht, exp(gg.ih),'r');
title('exponentiated post-spike kernel');
xlabel('time since spike (frames)');
ylabel('gain');
axis tight;

% Errors in STA and ML estimate
Estim_Error = [subspace(flts(:,1),flts(:,2)), subspace(flts(:,1),flts(:,3))] % In radians
