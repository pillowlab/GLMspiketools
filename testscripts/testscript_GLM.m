% testscript_GLM.m
%
% Test code for simulating and fitting the vanilla GLM: 
% - 1D (temporal) filter only, exponential nonlinearity
%
% Code Blocks:  
%   1. Set up model params
%   2. Show simulated model responses to stimuli
%   3. Make training data (for fitting params to data)
%   4. Fit GLM params via maximum-likelihood (requires optimization toolbox)

global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100; 

%% 1.  Set parameters and display for GLM % =============================

DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 20;    % Number of time bins in stimulus filter
ttk = [-nkt+1:0]';  % time relative to spike of stim filter taps
ggsim = makeSimStruct_GLM(nkt,DTsim); % Create GLM structure with default params

% === Make Fig: model params =======================
figure(1);  
subplot(2,2,1);
plot(ttk, ggsim.k);
title('stimulus kernel');
xlabel('time before spike (frames)');

subplot(2,2,2); % --------
plot(ggsim.iht, ggsim.iht*0, 'k--', ggsim.iht, ggsim.ih);
title('post-spike kernel h');
set(gca, 'ylim',[min(ggsim.ih)*1.1 max(ggsim.ih)*1.5]);

subplot(2,2,3); % --------
plot(ggsim.iht, exp(ggsim.iht*0), 'k--', ggsim.iht, exp(ggsim.ih));
title('exponentiated post-spike kernel h');
xlabel('time after spike (frames)');
ylabel('gain');

subplot(2,2,4); % --------
[iht,ihbasOrthog,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,DTsim);
plot(ggsim.iht, ihbasis);
axis tight;
xlabel('time after spike (frames)');
title('basis for h');


%% 2. Make GWN stimulus & simulate the glm model response. ===========

slen = 50; % Stimulus length (frames) 
swid = 1;  % Stimulus width  (pixels).  Must match # pixels in stim filter
Stim = randn(slen,swid);  % Gaussian white noise stimulus
[tsp, vmem,Ispk] = simGLM(ggsim, Stim);  % Simulate GLM response

% ==== Make Figure: repeat responses ========
figure(2); 
tt = [DTsim:DTsim:slen]';
subplot(221); %------------------------
plot(1:slen, Stim, 'k', 'linewidth', 2); 
title('GWN stimulus');
xlabel('time');
axis tight;

subplot(222); %------------------------
if ~isempty(tsp)
    plot(tt, vmem, tsp, max(vmem)*ones(size(tsp)), 'r.');
    title('net voltage (dots = spike times)');
    axis tight;
end

% -----Run a few repeat simulations ------
nrpts = 5;        % number of repeats to draw
subplot(223); %------------------------
plot(tt, vmem-Ispk, tt, Ispk, 'r');
title('stim (blue) and spike (red) -induced currents'); axis tight;
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
slen = 2500; % Stimulus length (frames);  More samples gives better fit
Stim = round(rand(slen,swid))*2-1;  %  Run model on long, binary stimulus
[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model
nsp = length(tsp);

% Compute STA and use as initial guess for k
sta0 = simpleSTC(Stim,tsp,nkt);
sta = reshape(sta0,nkt,[]);

% % -----------
% % Make param object with "true" params;
% ggTrue = makeFittingStruct_GLM(ggsim.k,DTsim,ggsim);
% ggTrue.tsp = tsp;
% ggTrue.tspi = 1;  % 1st spike to use for computing likelihood (eg, can ignore 1st n spikes)
% 
% % Check that conditional intensity calc is correct 
% % (if desired, compare rrT computed here with vmem computed above).
% [logliTrue, rrT,tt] = neglogli_GLM(ggTrue,Stim);
% % -----------

%% 4. Do ML fitting of params with simulated data %=====================

%  Initialize params for fitting --------------
gg0 = makeFittingStruct_GLM(sta,DTsim);  % projects sta into basis for fitting k
gg0.tsp = tsp;  % Insert spikes into fitting struct
gg0.tspi = 1;   % First spike to use (you can ask it to ignore the first "n" spikes)
[logli0,rr0,tt] = neglogli_GLM(gg0,Stim); % Compute logli of initial params (if desired)

% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 100};
[gg, negloglival] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% 5. Plot results ======================================================
figure(3);
ttk = -nkt+1:0;
subplot(221);  % True filter  % ---------------
plot(ttk, ggsim.k, 'k', ttk, sta, ttk, gg.k, 'r');
title('Stim filters (True=blck, STA=blue, ML=red)');

subplot(223);
flts = normalizecols([ggsim.k, sta, gg.k]);
plot(ttk, flts(:,1),'k', ttk,flts(:,2), ttk, flts(:,3), 'r');
title('Normalized Stim filters');
xlabel('time before spike (frames)');

subplot(222); % ----------------------------------
plot(ggsim.iht, ggsim.ih, gg.iht, gg.ihbas*gg.ih,'r');
title('post-spike kernel');
axis tight;

subplot(224); % ----------------------------------
plot(ggsim.iht, exp(ggsim.ih), gg.iht, exp(gg.ihbas*gg.ih),'r');
title('exponentiated post-spike kernel');
xlabel('time since spike (frames)');
ylabel('gain');
axis tight;

% Errors in STA and ML estimate
Estim_Error = [subspace(flts(:,1),flts(:,2)), subspace(flts(:,1),flts(:,3))] % In radians
