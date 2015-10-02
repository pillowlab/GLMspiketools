% testscript_GLM_coupled.m
%
% Test code for simulating and fitting a coupled GLM (2 neurons).
%
% Notes:
%   Fitting code uses same functions as for single-cell responses.
%   Simulation code requires new structures / functions
%     (due to the need to pass activity between neurons)
%
% Code Blocks:  
%   1. Set up model params and plot
%   2. Show samples from simulated model
%   3. Generate training data
%   4. Fit simulated dataset w/ maximum-likelihood (ML)
%   5. Plot results

global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100;


%%  1.  Set parameters and display for GLM  ============ %

DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 20;    % Number of time bins in filter;
ttk = [-nkt+1:0]'; 
ggsim1 = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
ggsim2 = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
k = ggsim1.k;  % Stimulus filter
ggsim2.k = shift(k,-3)*1.2;
ggsim = makeSimStruct_GLMcpl(ggsim1,ggsim2);

% Make some coupling kernels
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,DTsim);
hhcpl = ihbasis*[.6;.47;.25;0;0]*2;
hhcpl(:,2) = ihbasis*[-1;-1;0;0;.25]*2;
ggsim.ih(:,2,1) = hhcpl(:,2); % 2nd cell coupling to first
ggsim.ih(:,1,2) = hhcpl(:,1); % 1st cell coupling to second

% === Make Fig: model params =======================
figure(1);  
subplot(2,2,1);
plot(ttk, ggsim.k(:,:,1));
title('neuron 1 stim kernel');
subplot(2,2,3);
plot(ttk, ggsim.k(:,:,2));
title('neuron 2 stim kernel');
xlabel('time (frames)');

subplot(2,2,2); % --------
plot(ggsim.iht, exp(ggsim.ih(:,:,1)), ggsim.iht, ggsim.iht*0+1, 'k--');
title('spike kernels into Cell 1 (exponentiated)');
legend('from 1', 'from 2', 'location', 'northeast');
ylabel('gain');

subplot(2,2,4); % --------
plot(ggsim.iht, exp(ggsim.ih(:,:,2)), ggsim.iht, ggsim.iht*0+1, 'k--');
title('spike kernels into Cell 2 (exponentiated)');
legend('from 1', 'from 2', 'location', 'northeast');
ylabel('multiplicative gain');
xlabel('time after spike (frames)');


%% 2. Make GWN stimulus & simulate the glm model response. ========= %
% 
slen = 50; % Stimulus length (frames) & width (# pixels)
swid = size(ggsim.k,2);
Stim = 2*randn(slen,swid);  % Gaussian white noise stimulus
[tsp, vmem,Ispk] = simGLM(ggsim, Stim); % Simulate GLM response

% ==== Make Figure ========
figure(2); 
tt = [DTsim:DTsim:slen]';
subplot(321); %------------------------
plot(1:slen, Stim, 'k', 'linewidth', 2); 
title('GWN stimulus');
axis tight;

subplot(3,2,3); %------------------------
plot(tt, vmem(:,1), tsp{1}, max(vmem(:,1))*ones(size(tsp{1})), 'ro');
title('cell 1: net voltage + spikes');
axis tight;
% -----Run repeat simulations ------
nrpts = 5;        % number of repeats to draw
ylabel('filter output');

subplot(3,2,4); %------------------------
plot(tt, vmem(:,2), 'b', tsp{2}, max(vmem(:,2))*ones(size(tsp{2})), 'ro');
title('cell 2: net voltage + spikes');
axis tight;

subplot(325)
plot(tt, vmem(:,1)-Ispk(:,1), 'k', tt, Ispk(:,1), 'r');
title('stim-induced & spike-induced currents'); 
axis tight;
xlabel('time (frames)');
ylabel('filter output');

subplot(326)
plot(tt, vmem(:,2)-Ispk(:,2), 'k', tt, Ispk(:,2), 'r');
title('stim-induced & spike-induced currents'); 
axis tight;
xlabel('time (frames)');


%% 3. Generate some training data
slen = 2500;  % Stimulus length (frames);  More samples gives better fit
Stim = round(rand(slen,swid))*4-2;  %  Run model on long, binary stimulus
[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model

% -------------- Compute STAs------------
nsp = length(tsp{1});
sta0 = simpleSTC(Stim,tsp{1},nkt); % Compute STA 1
sta1 = reshape(sta0,nkt,[]); 

sta0 = simpleSTC(Stim,tsp{2},nkt); % Compute STA 2
sta2 = reshape(sta0,nkt,[]); 


% % ---------------
% % Make param object with "true" params;
% ggTrue1 = makeFittingStruct_GLM(ggsim.k(:,:,1),DTsim,ggsim,1);
% ggTrue1.tsp = tsp{1}; % cell 1 spike times (vector) 
% ggTrue1.tspi = 1; % 1st spike to use for computing likelihood (eg, can ignore 1st n spikes)
% ggTrue1.tsp2 = tsp(2); % spike trains from "coupled" cells (cell array of vectors)
% 
% % Check that conditional intensity calc is correct:
% fprintf('Fitting first neuron\n');
% [logliTrue, rrT,tt] = neglogli_GLM(ggTrue1,Stim);
% % ---------------


%% 4. Do ML fitting of GLM params: first cell ============= %
%  Initialize params for fitting --------------
gg0 = makeFittingStruct_GLM(sta1,DTsim,ggsim,1);  % Initialize params for fitting struct w/ sta
gg0.ih = gg0.ih*0;  % Initialize to zero
gg0.dc = gg0.dc*0;  % Initialize to zero

gg0.tsp = tsp{1};   % cell 2 spike times (vector)
gg0.tsp2 = tsp(2);  % spike trains from "coupled" cells (cell array of vectors)
gg0.tspi = 1; % 1st spike to use for computing likelihood (eg, can ignore 1st n spikes)

% Do ML estimation of model params
fprintf('Fitting first neuron\n');
opts = {'display', 'iter', 'maxiter', 100};
[gg1, negloglival1] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% Fit second cell
gg0b = gg0; % initial parameters for fitting 
gg0b.tsp = tsp{2};   % cell 2 spike times (vector)
gg0b.tsp2 = tsp(1);  % spike trains from "coupled" cells (cell array of vectors)
gg0b.kt = inv(gg0.ktbas'*gg0.ktbas)*gg0.ktbas'*sta2; % Project STA2 into basis 
gg0b.k = gg0b.ktbas*gg0b.kt; % Project STA onto basis for fitting
fprintf('\nFitting second neuron\n');
[gg2, negloglival2] = MLfit_GLM(gg0b,Stim,opts); % do ML (requires optimization toolbox)


%% --- Plot results ----------------------------
figure(3);
ttk = -nkt+1:0;
subplot(221);  % Filters cell 1 % ---------------
plot(ttk, ggsim1.k, 'k', ttk, gg1.k, 'r');
title('Cell 1 stim filt (True=blck, ML=red)');

subplot(223);  % Filters cell 2 % ---------------
plot(ttk, ggsim2.k, 'k', ttk, gg2.k, 'r');
title('Cell 2 stim filt (True=blck, ML=red)');
xlabel('time (frames)')

subplot(222); % ----------------------------------
plot(ggsim.iht, exp(ggsim.ih(:,:,1)), 'k', gg1.iht, exp(gg1.ihbas*[gg1.ih gg1.ih2]));
title('Cell 1: exponentiated post-spk kernels');
axis tight;

subplot(224); % ----------------------------------
plot(ggsim.iht, exp(ggsim.ih(:,:,2)), 'k', gg2.iht, exp(gg2.ihbas*[gg2.ih gg2.ih2]));
title('Cell 2: exponentiated post-spk kernels');
xlabel('time (frames)')
axis tight;

