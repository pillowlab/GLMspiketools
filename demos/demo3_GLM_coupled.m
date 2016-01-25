% testscript3_GLM_coupled.m
%
% Demo script for simulating and fitting a coupled GLM (2 neurons).
%
% Notes:
%  - Fitting code uses same functions as for single-cell responses.
%  -  Simulation code requires new structures / functions
%     (due to the need to pass activity between neurons)

% Make sure paths are set (assumes this script called from 'demos' directory)
cd ..; setpaths; cd demos/


%%  ===== 1. Set parameters and display for GLM  ============ %

dtStim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
dtSp = .001;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
nkt = 20;    % Number of time bins in stimulus filter k
ggsim1 = makeSimStruct_GLM(nkt,dtStim,dtSp);  % Create GLM struct with default params
ggsim2 = makeSimStruct_GLM(nkt,dtStim,dtSp);  % Create GLM struct with default params

% Change second neuron's stimulus filter
[ktbas,ktbasis] = makeBasis_StimKernel(ggsim2.ktbasprs,nkt);
ggsim2.k = ktbasis*[0 0 .1 .25 .5 .5 .25 -.25 -.5 -.25]'; % more delayed filter
ggsim = makeSimStruct_GLMcpl(ggsim1,ggsim2);

% Make some coupling kernels
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,dtSp);
hhcpl = ihbasis*[.25;.25;.125;0;0];
hhcpl(:,2) = ihbasis*[-2;-1;0;.25;.25];
ggsim.ih(:,2,1) = hhcpl(:,2); % 2nd cell coupling to first
ggsim.ih(:,1,2) = hhcpl(:,1); % 1st cell coupling to second

% === Make Fig: model params =======================
ttk = (-nkt+1:0)'*dtStim; 
subplot(2,2,1); % -------------------------------
plot(ttk, ggsim.k(:,:,1));
title('neuron 1 stim kernel');
subplot(2,2,3); % -------------------------------
plot(ttk, ggsim.k(:,:,2));
title('neuron 2 stim kernel'); xlabel('time (frames)');
subplot(2,2,2); % --------------------------------
plot(ggsim.iht, exp(ggsim.ih(:,:,1)), ggsim.iht, ggsim.iht*0+1, 'k--');
title('spike kernels into Cell 1 (exponentiated)');
legend('from 1', 'from 2', 'location', 'northeast'); ylabel('gain');
subplot(2,2,4); % ---------------------------------
plot(ggsim.iht, exp(ggsim.ih(:,:,2)), ggsim.iht, ggsim.iht*0+1, 'k--');
title('spike kernels into Cell 2 (exponentiated)');
legend('from 1', 'from 2', 'location', 'northeast');
ylabel('gain'); xlabel('time after spike (frames)');


%% ===== 2. Run short simulation for visualization purposes ========= %
% 
slen = 50; % Stimulus length (frames) & width (# pixels)
swid = size(ggsim.k,2); % stimulus width
Stim = 2*randn(slen,swid);  % Gaussian white noise stimulus
[tsp,~,Itot,Istm] = simGLM(ggsim, Stim); % Simulate GLM response
Isp = Itot-Istm; % net spike-history output

% ==== Plot some traces of simulated response  ========
tt = (dtSp:dtSp:slen*dtStim)';
subplot(321); %------------------------
plot(1:slen, Stim, 'k', 'linewidth', 2); 
title('GWN stimulus'); axis tight;
subplot(3,2,3); %------------------------
plot(tt, Itot(:,1),'b',tsp{1}, max(Itot(:,1))*ones(size(tsp{1})), 'ro');
title('cell 1: net voltage + spikes');
axis([0 slen*dtStim, min(Itot(:,1)) max(Itot(:,1))*1.01]); ylabel('filter output');
subplot(3,2,4); %------------------------
plot(tt, Itot(:,2), 'b', tsp{2}, max(Itot(:,2))*ones(size(tsp{2})), 'ro');
title('cell 2: net voltage + spikes');
axis([0 slen*dtStim, min(Itot(:,2)) max(Itot(:,2))*1.01]); 
subplot(325)  % --------------------------
plot(tt, Istm(:,1), 'k', tt, Isp(:,1), 'r');
title('stim-induced & spike-induced currents');  axis tight;
xlabel('time (frames)'); ylabel('filter output');
subplot(326); %---------------------------
plot(tt, Istm(:,2), 'k', tt, Isp(:,2), 'r');
title('stim-induced & spike-induced currents'); 
axis tight; xlabel('time (frames)');


%% ===== 3. Generate some training data =============================== %%

slen = 20000;  % Stimulus length (frames);  More samples gives better fit
Stim = round(rand(slen,swid))-.5;  %  Run model on long, binary stimulus
[tsp,sps,Itot,ispk] = simGLM(ggsim,Stim);  % run model


%% ===== 4. Fit cell #1 (with coupling from cell #2) =================== %%

% Compute Spike-triggered averages
sps1 = sum(reshape(sps(:,1),[],slen),1)'; % bin spikes in bins the size of stimulus
sta1 = simpleSTC(Stim,sps1,nkt); % Compute STA 1
sta1 = reshape(sta1,nkt,[]); 

sps2 = sum(reshape(sps(:,2),[],slen),1)'; % rebinned spike train
sta2 = simpleSTC(Stim,sps2,nkt); % Compute STA 2
sta2 = reshape(sta2,nkt,[]); 

% Initialize param struct for fitting 
gg0 = makeFittingStruct_GLM(dtStim,dtSp);  % Initialize params for fitting struct 

% Initialize fields (using h and k bases computed above)
gg0.ktbas = ktbas; % k basis
gg0.ihbas = ihbas; % h self-coupling basis
gg0.ihbas2 = ihbas; % h coupling-filter basis
nktbasis = size(ktbas,2); % number of basis vectors in k basis
nhbasis = size(ihbas,2); % number of basis vectors in h basis
gg0.kt = 0.1*(ktbas\sta1); % initial params from scaled-down sta 
gg0.k = gg0.ktbas*gg0.kt;  % initial setting of k filter
gg0.ihw = zeros(nhbasis,1); % params for self-coupling filter
gg0.ihw2 = zeros(nhbasis,1); % params for cross-coupling filter
gg0.ih = [gg0.ihbas*gg0.ihw gg0.ihbas2*gg0.ihw2];
gg0.iht = iht;
gg0.dc = 0; % Initialize dc term to zero
gg0.couplednums = 2; % number of cell coupled to this one (for clarity)

% Set spike responses for cell 1 and coupled cell
gg0.sps = sps(:,1);  
gg0.sps2 = sps(:,2); 

% Compute initial value of negative log-likelihood (just to inspect)
[neglogli0,rr] = neglogli_GLM(gg0,Stim);

% Do ML fitting
fprintf('Fitting neuron 1:  initial neglogli0 = %.3f\n', neglogli0);
opts = {'display', 'iter', 'maxiter', 100};
[gg1, neglogli1] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% ===== 5. Fit cell #2 (with coupling from cell #1) ==================

gg0b = gg0; % initial parameters for fitting 
gg0b.sps = sps(:,2); % cell 2 spikes
gg0b.sps2 = sps(:,1); % spike trains from coupled cells 
gg0.kt = 0.1*(ktbas\sta2); % initial params from scaled-down sta 
gg0b.k = gg0b.ktbas*gg0b.kt; % Project STA onto basis for fitting
gg0.couplednums = 1; % number of cell coupled to this one

% Compute initial value of negative log-likelihood (just to inspect)
[neglogli0b] = neglogli_GLM(gg0b,Stim); % initial value of negative logli

% Do ML fitting
fprintf('Fitting neuron 2: initial neglogli = %.3f\n', neglogli0b);
[gg2, neglogli2] = MLfit_GLM(gg0b,Stim,opts); % do ML (requires optimization toolbox)


%% ===== 6. Plot fits  ============================================= %%

ttk = (-nkt+1:0)*dtStim; % time indices for k

subplot(221);  % Stim filters cell 1 % ---------------
plot(ttk, ggsim1.k(:,1), 'k', ttk, gg1.k, 'r');
title('Cell 1: k filter');
legend('true', 'estim', 'location', 'northwest');
subplot(223);  % Stim filters cell 2 % ---------------
plot(ttk, ggsim2.k, 'k', ttk, gg2.k, 'r');
title('Cell 2: k filter');
xlabel('time (s)')
subplot(222); % --Spike filters cell 1 % -------------
plot(ggsim.iht, (ggsim.ih(:,:,1)), gg1.iht, (gg1.ih), '--');
title('exponentiated incoming h filters');
legend('true h11', 'true h21', 'estim h11', 'estim h21');
axis tight;
subplot(224); % --Spike filters cell 2 % ------------- 
plot(ggsim.iht, (ggsim.ih(:,:,2)), gg2.iht, (gg2.ih), '--');
title('exponentiated incoming h filters');
axis tight; xlabel('time (s)')
legend('true h22', 'true h12', 'estim h22', 'estim h12');


