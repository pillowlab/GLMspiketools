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


%% 1.  Set parameters and display for GLM % =============================

dtStim = .01; % Bin size for stimulus (in seconds).  (Equiv to 100Hz frame rate)
dtSp = .001;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
nkt = 30;    % Number of time bins in stimulus filter k
ttk = dtStim*(-nkt+1:0)';  % time relative to spike of stim filter taps
ggsim = makeSimStruct_GLM(nkt,dtStim,dtSp); % Create GLM structure with default params

% === Make Fig: model params =======================
%figure(1);  
subplot(2,2,1);
plot(ttk, ggsim.k);
title('stimulus kernel');
xlabel('time (s)');

subplot(2,2,2); % --------
[iht,ihbasOrthog,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,dtSp);
plot(ggsim.iht, ihbasis); 
title('basis for h'); axis tight;

subplot(2,2,3); % --------
plot(ggsim.iht, exp(ggsim.iht*0), 'k--', ggsim.iht, exp(ggsim.ih));
title('exponentiated h');
xlabel('time after spike (s)');
ylabel('gain'); axis tight;

subplot(2,2,4); % --------
plot(ggsim.iht, ggsim.iht*0, 'k--', ggsim.iht, ggsim.ih);
title('post-spike kernel h'); axis tight;
xlabel('time after spike (s)');
set(gca, 'ylim',[min(ggsim.ih)*1.1 max(ggsim.ih)*1.5]);




%% 2. Make GWN stimulus & simulate the glm model response. ===========

slen = 50; % Stimulus length (frames) 
swid = 1;  % Stimulus width  (pixels).  Must match # pixels in stim filter
Stim = randn(slen,swid);  % Gaussian white noise stimulus
[tsp, vmem,Ispk] = simGLM(ggsim, Stim);  % Simulate GLM response

% ==== Make Figure: repeat responses ========
%figure(2); 
tt = (dtSp:dtSp:(slen*dtStim))';
subplot(221); %------------------------
plot(1:slen, Stim, 'k', 'linewidth', 2); 
title('GWN stimulus');
xlabel('time');
axis tight;

subplot(222); %------------------------
if ~isempty(tsp)
    plot(tt, vmem, tsp, max(vmem)*ones(size(tsp)), 'r.');
    title('net voltage (dots = spike times)');
    axis([0 tt(end), min(vmem) max(vmem*1.01)]);
end

% -----Run a few repeat simulations ------
nrpts = 5;        % number of repeats to draw
subplot(223); %------------------------
plot(tt, vmem-Ispk, tt, Ispk, 'r');
title('stim (blue) and spike (red) -induced currents'); 
axis tight;
xlabel('time (frames)');
subplot(224);
for j = 1:nrpts;
    [tsp1,vmem1] = simGLM(ggsim,Stim);
    plot(tt,vmem1, 'color', rand(3,1)); hold on;
end
axis tight; hold off;
title('repeated responses to same stim');
xlabel('time (frames)');


%% 3. Make some training data  %========================================
slen = 5000; % Stimulus length (frames);  More samples gives better fit
Stim = round(rand(slen,swid))*2-1;  %  Run model on long, binary stimulus

[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model
nsp = length(tsp);

% Compute STA and use as initial guess for k
sta0 = simpleSTC(Stim,tsp,nkt);
sta = reshape(sta0,nkt,[]);

exptmask= [100 slen];  % data range to use for fitting

% -----------
% Make param object with "true" params
ggTrue = makeFittingStruct_GLM(ggsim.k,dtStim,dtSp,ggsim);
ggTrue.tsp = tsp;
ggTrue.mask = exptmask;

% Check that conditional intensity calc is correct 
% (if desired, compare rrT computed here with vmem computed above).
[logliTrue, rrT,tt] = neglogli_GLM(ggTrue,Stim);
subplot(211); plot(tt,vmem,tt,log(rrT));
title('total filter output (computed 2 ways)');
subplot(212); plot(tt,log(rrT)-vmem);
title('difference');
% -----------

%% 4. Do ML fitting of params with simulated data %=====================

%  Initialize params for fitting --------------
gg0 = makeFittingStruct_GLM(sta,dtSp);  % projects sta into basis for fitting k
gg0.tsp = tsp;  % Insert spikes into fitting struct
gg0.mask = exptmask;
[logli0,rr0,tt] = neglogli_GLM(gg0,Stim); % Compute logli of initial params (if desired)
fprintf('Initial value of negative log-li: %.3f\n', logli0);

% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 100}; 
[gg, negloglival] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)

%% 5. Plot results ======================================================
%figure(3);
ttk = -nkt+1:0;
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
