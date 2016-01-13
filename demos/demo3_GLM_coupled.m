% testscript3_GLM_coupled.m
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
%   4. Fit simulated dataset via maximum-likelihood (ML)
%   5. Plot results

global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
RefreshRate = 100;


%%  1.  Set parameters and display for GLM  ============ %

DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
nkt = 20;    % Number of time bins in filter;
ttk = [-nkt+1:0]'; 
ggsim1 = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
ggsim2 = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params

% Change second neuron's stimulus filter
[ktbas,ktbasis] = makeBasis_StimKernel(ggsim2.ktbasprs,nkt);
ggsim2.k = ktbasis*[0 0 .1 .25 .5 .5 .25 -.25 -.5 -.25]'; % more delayed filter
ggsim = makeSimStruct_GLMcpl(ggsim1,ggsim2);

% Make some coupling kernels
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,DTsim);
hhcpl = ihbasis*[.5;.5;.25;0;0];
hhcpl(:,2) = ihbasis*[-1;-1;0;.25;.25];
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


%% 2. Make Gaussian White Noise stimulus & simulate glm resp ========= %
% 
slen = 50; % Stimulus length (frames) & width (# pixels)
swid = size(ggsim.k,2); % stimulus width
Stim = 2*randn(slen,swid);  % Gaussian white noise stimulus
[tsp, vmem,Ispk] = simGLM(ggsim, Stim); % Simulate GLM response

% ==== Plot some traces of simulated response  ========
figure(2); 
tt = [DTsim:DTsim:slen]';
subplot(321); %------------------------
plot(1:slen, Stim, 'k', 'linewidth', 2); 
title('GWN stimulus');
axis tight;

subplot(3,2,3); %------------------------
plot(tt, vmem(:,1), tsp{1}, max(vmem(:,1))*ones(size(tsp{1})), 'ro');
title('cell 1: net voltage + spikes');
axis([0 slen, min(vmem(:,1)) max(vmem(:,1))*1.01]);
% -----Run repeat simulations ------
nrpts = 5;        % number of repeats to draw
ylabel('filter output');

subplot(3,2,4); %------------------------
plot(tt, vmem(:,2), 'b', tsp{2}, max(vmem(:,2))*ones(size(tsp{2})), 'ro');
title('cell 2: net voltage + spikes');
axis([0 slen, min(vmem(:,2)) max(vmem(:,2))*1.01]);
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

exptmask= [100 slen];  % data range to use for fitting

% %% ---------------
% % Uncomment to check consistency in computing conditional intensity between
% % simulation and fitting code
% % ----------------
% 
% % % Make param object with "true" params;  
% cellnum = 1;
% ggTrue1 = makeFittingStruct_GLM(ggsim.k(:,:,1),DTsim,ggsim,cellnum);
% ggTrue1.tsp = tsp{1}; % cell 1 spike times (vector) 
% ggTrue1.tsp2 = tsp(2); % spike trains from "coupled" cells (cell array of spike-time vectors)
% ggTrue1.mask = exptmask;
% % Check that conditional intensity calc is correct:
% [logliTrue1, rrT1,tt] = neglogli_GLM(ggTrue1,Stim);
% subplot(211); plot(tt,vmem(:,1),tt,log(rrT1));
% title('total linear filter output');
% subplot(212); plot(tt,log(rrT1)-vmem(:,1));
% title('difference');
% 
% %%
% cellnum = 2;
% ggTrue2 = makeFittingStruct_GLM(ggsim.k(:,:,2),DTsim,ggsim,cellnum);
% ggTrue2.tsp = tsp{2};  % cell 2 spike times (vector) 
% ggTrue2.tsp2 = tsp(1); % spike trains from "coupled" cells (cell array of spike-time vectors)
% ggTrue2.mask = exptmask;
% % Check that conditional intensity calc is correct:
% [logliTrue2, rrT2,tt] = neglogli_GLM(ggTrue2,Stim);
% subplot(211); plot(tt,vmem(:,2),tt,log(rrT2));
% title('total linear filter output');
% subplot(212); plot(tt,log(rrT2)-vmem(:,2));
% title('difference');
% % % ---------------


%% 4. Do ML fitting of params with simulated data ============= %
%  Initialize params for fitting --------------
gg0 = makeFittingStruct_GLM(sta1,DTsim,ggsim,1);  % Initialize params for fitting struct w/ sta
gg0.ihw = gg0.ihw*0;  % Initialize to zero
gg0.ihw2 = gg0.ihw2*0;  % Initialize to zero
gg0.dc = gg0.dc*0;  % Initialize to zero

gg0.tsp = tsp{1};   % cell 2 spike times (vector)
gg0.tsp2 = tsp(2);  % spike trains from "coupled" cells (cell array of vectors)
gg0.mask = exptmask;

[logli0] = neglogli_GLM(gg0,Stim);

% Do ML estimation of model params
fprintf('Fitting first neuron ( logli0=%.3f )\n', logli0);
opts = {'display', 'iter', 'maxiter', 100};
[gg1, negloglival1] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)


%% Fit other cell
gg0b = gg0; % initial parameters for fitting 
gg0b.tsp = tsp{2};   % cell 2 spike times (vector)
gg0b.tsp2 = tsp(1);  % spike trains from "coupled" cells (cell array of vectors)
gg0b.kt = inv(gg0.ktbas'*gg0.ktbas)*gg0.ktbas'*sta2; % Project STA2 into basis 
gg0b.k = gg0b.ktbas*gg0b.kt; % Project STA onto basis for fitting

[logli0b] = neglogli_GLM(gg0b,Stim);
fprintf('Fitting second neuron (logli0=%.3f)\n', logli0b);
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
plot(ggsim.iht, exp(ggsim.ih(:,:,1)), 'k', gg1.iht, exp(gg1.ih));
title('Cell 1: exponentiated post-spk kernels');
axis tight;

subplot(224); % ----------------------------------
plot(ggsim.iht, exp(ggsim.ih(:,:,2)), 'k', gg2.iht, exp(gg2.ih));
title('Cell 2: exponentiated post-spk kernels');
xlabel('time (frames)')
axis tight;

