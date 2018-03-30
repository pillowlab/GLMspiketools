% testscript3_GLM_coupled.m
%
% Demo script for simulating and fitting a coupled GLM (2 neurons).
%
% Notes:
%  - Fitting code uses same functions as for single-cell responses.
%  -  Simulation code requires new structures / functions
%     (due to the need to pass activity between neurons)

% Make sure paths are set (assumes this script called from 'demos' directory)
cd ..; setpaths_GLMspiketools; cd demos/


%%  ===== 1. Set parameters for simulating a GLM  ============ %

dtStim = .001; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
dtSp = .001;  % Bin size for simulating model & computing likelihood (must evenly divide dtStim);
nkt = 1;    % Number of time bins in stimulus filter k %% FIX THIS SO IT DOESN'T MAKE NAN
gg = makeSimStruct_GLM(nkt,dtStim,dtSp);  % Create GLM struct with default params

%% Set up basis for self-coupling filters

% Make basis for self-coupling term
ihbasprs.ncols = 3; % number of basis vectors
ihbasprs.hpeaks = [.002, .005]; % peak of 1st and last basis vector
ihbasprs.b = .001;  % scaling (smaller -> more logarithmic scaling)
ihbasprs.absref = .002; % absolute refractory period basis vector (optional)
% Make basis 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dtSp);
nht = length(iht); % number of bins

% Make basis for cross-coupling term
ihbasprs2.ncols = 1;  % number of basis vectors
ihbasprs2.hpeaks = [0.001,.005]; % put peak at 10ms and "effective" 1st peak at 0
ihbasprs2.b = .002;  % smaller -> more logarithmic scaling
ihbasprs2.absref = []; % no abs-refracotry period for this one
% Make basis
[iht2,ihbas2,ihbasis2] = makeBasis_PostSpike(ihbasprs2,dtSp);
nht2 = length(iht2);

% pad to put them into the same time bins, if necessary
if nht2>nht
    % padd ih1 with zeros
    iht = iht2; zz = zeros(nht2-nht,ihbasprs.ncols);
    ihbas = [ihbas;zz]; ihbasis = [ihbasis;zz]; nht=nht2;
elseif nht2<nht
    % padd ih1 with zeros
    iht2 = iht; zz = zeros(nht-nht2,ihbasprs2.ncols);
    ihbas2 = [ihbas2;zz]; ihbasis2 = [ihbasis2;zz]; nht2=nht;
end    

% plot it
subplot(211); 
plot(iht, ihbasis, ihbasprs.hpeaks, 1, '*', 'linewidth', 2);
xlabel('time after spike (ms)'); title('post-spike basis');
subplot(212); 
plot(iht2, ihbasis2, ihbasprs2.hpeaks, 1, '*', 'linewidth', 2);
xlabel('time after spike (ms)'); title('coupling basis');


%% Set self-coupling weights

wself = [-5; .2; -.15]; % weights for self-coupling term
ihself = ihbasis*wself; % self-coupling filter
wcpl = 0.5; % weights for cross-coupling term
ihcpl = ihbasis2*wcpl; % cross-coupling filter
clf; plot(iht, exp(ihself), iht, exp(ihcpl), iht, iht*0+1, 'k--');
legend('self-coupling', 'cross-coupling');
xlabel('time lag (s)');
ylabel('gain (sp/s)');


%% Set up multi-neuron GLM

nneur = 3;
k = [.85 .9 0.7]; % stimulus weights
gg.k = permute(k,[1,3,2]);  % stimulus weights

gg.iht = iht;
gg.ih = zeros(nht,nneur,nneur);
gg.ih(:,:,1) = [ihself, 1.25*ihcpl, -1*ihcpl]; % input weights to neuron 1
gg.ih(:,:,2) = [1.8*ihcpl, 1.5*ihself, -.5*ihcpl]; % input weights to neuron 2
gg.ih(:,:,3) = [1.5*ihcpl, .5*ihcpl, ihself]; % input weights to neuron 3


%% ===== 2. Run short simulation for visualization purposes ========= %

slen = 50; % Stimulus length (frames) & width (# pixels)
swid = 1; % width of stimulus
stimsd = 1;  % contrast of stimulus
Stim = stimsd*randn(slen,swid);  % Gaussian white noise stimulus
[tsp,~,Itot,Istm] = simGLM(gg, Stim); % Simulate GLM response
Isp = Itot-Istm; % net spike-history output

% ==== Plot some traces of simulated response  ========
tt = (dtSp:dtSp:slen*dtStim)';
subplot(131)  % % ==== neuron 1 ===========
plot(tt, Istm(:,1), 'k', tt, Isp(:,1), 'r',  tt, tt*0, 'k--', ...
    tsp{1}, ones(size(tsp{1})), 'ro');
title('cell 1');  axis tight;
xlabel('time (s)'); ylabel('filter output');
legend('stim filter + dc', 'spike-hist filter');
subplot(132)  % ==== neuron 2 ===========
plot(tt, Istm(:,2), 'k', tt, Isp(:,2), 'r',  tt, tt*0, 'k--', ...
    tsp{2}, ones(size(tsp{2})), 'ro');
title('cell 2');  axis tight;
xlabel('time (s)'); 
subplot(133)  % ==== neuron 3 ===========
plot(tt, Istm(:,3), 'k', tt, Isp(:,3), 'r',  tt, tt*0, 'k--', ...
    tsp{3}, ones(size(tsp{3})), 'ro');
title('cell 3');  axis tight;
xlabel('time (s)'); 

%% ===== 3. Generate some training data =============================== %%

slen = 2e5;  % Stimulus length (in bins);  
Stim = stimsd*randn(slen,swid);  % Gaussian white noise stimulus
[tsp,sps,Itot,ispk] = simGLM(gg,Stim);  % run model


%% ===== 4. Fit cell #1 (with coupling from cell #2 and #3) =================== %%


% Initialize param struct for fitting 
gg1in = makeFittingStruct_GLM(dtStim,dtSp);  % Initialize params for fitting struct 

% Initialize fields (using h bases computed above)
gg1in.ktbas = 1; % k basis
gg1in.ihbas = ihbas; % h self-coupling basis
gg1in.ihbas2 = ihbas2; % h coupling-filter basis
nktbasis = 1; % number of basis vectors in k basis
nhbasis = size(ihbas,2); % number of basis vectors in h basis
nhbasis2 = size(ihbas2,2); % number of basis vectors in h basis
gg1in.kt = 1; % initial params from scaled-down sta 
gg1in.k = 1;  % initial setting of k filter
gg1in.ihw = zeros(nhbasis,1); % init params for self-coupling filter
gg1in.ihw2 = zeros(nhbasis2,nneur-1); % init params for cross-coupling filter
gg1in.ih = [gg1in.ihbas*gg1in.ihw gg1in.ihbas2*gg1in.ihw2];
gg1in.iht = iht;
gg1in.dc = 0; % Initialize dc term to zero

% Set fields for fitting cell #1
couplednums = [2 3];  % the cells coupled to this one
gg1in.couplednums = couplednums; % cell numbers of cells coupled to this one 
gg1in.sps = sps(:,1);  % Set spike responses for cell 1 
gg1in.sps2 = sps(:,couplednums); % spikes from coupled cells

% Compute initial value of negative log-likelihood (just to inspect)
[neglogli0,rr] = neglogli_GLM(gg1in,Stim);

% Do ML fitting
fprintf('Fitting neuron 1:  initial neglogli0 = %.3f\n', neglogli0);
opts = {'display', 'iter', 'maxiter', 100};
[gg1, neglogli1] = MLfit_GLM(gg1in,Stim,opts); % do ML (requires optimization toolbox)

%% ===== 5. Fit cell #2 (with coupling from cell #1 and #3) ==================

cellnum = 2;
couplednums = setdiff(1:3, cellnum); % the cells coupled to this one

gg2in = gg1in; % initial parameters for fitting 
gg2in.sps = sps(:,cellnum); % cell 2 spikes
gg2in.sps2 = sps(:,couplednums); % spike trains from coupled cells 
gg2in.couplednums = couplednums; % cells coupled to this one

% Do ML fitting
fprintf('Fitting neuron 2\n');
[gg2, neglogli2] = MLfit_GLM(gg2in,Stim,opts); % do ML (requires optimization toolbox)


%% ===== 6. Fit cell #3 (with coupling from cell #1 and #2) ==================

cellnum = 3;
couplednums = setdiff(1:3, cellnum);

gg3in = gg1in; % initial parameters for fitting 
gg3in.sps = sps(:,cellnum); % cell 3 spikes
gg3in.sps2 = sps(:,couplednums); % spike trains from coupled cells 
gg3in.couplednums = couplednums; % cells coupled to this one

% Do ML fitting
fprintf('Fitting neuron 3\n');
[gg3, neglogli3] = MLfit_GLM(gg3in,Stim,opts); % do ML (requires optimization toolbox)


%% ===== 6. Plot fits  ============================================= %%

colors = get(gca,'colororder');
set(gcf,'defaultAxesColorOrder',colors(1:3,:)); % use only 3 colors
lw = 2; % linewidth
ymax = max(exp([gg1.ih(:);gg2.ih(:);gg3.ih(:);gg.ih(:)])); % max of y range

subplot(131); % --Spike filters cell 1 % -------------
plot(gg.iht, exp(gg.ih(:,:,1)), gg1.iht, exp((gg1.ih)), '--', 'linewidth', lw);
hold on; plot(gg.iht, gg.iht*0+1, 'k'); hold off;
legend('true h11', 'true h21', 'true h31', 'estim h11', 'estim h21', 'estim h31', 'location', 'southeast');
title('incoming filters: cell 1'); axis tight; set(gca,'ylim',[0,ymax]); 
ylabel('gain (sp/s)'); xlabel('time after spike (s)');

subplot(132); % --Spike filters cell 2 % -------------
plot(gg.iht, exp(gg.ih(:,[2 1 3],2)), gg2.iht, exp(gg2.ih), '--', 'linewidth', lw);
hold on; plot(gg.iht, gg.iht*0+1, 'k'); hold off;
legend('true h22', 'true h12', 'true h32', 'estim h22', 'estim h12', 'estim h32', 'location', 'southeast');
title('cell 2'); axis tight; set(gca,'ylim',[0,ymax]);

subplot(133); % --Spike filters cell 3 % -------------
plot(gg.iht, exp(gg.ih(:,[3 1 2],3)), gg3.iht, exp(gg3.ih), '--', 'linewidth', lw);
hold on; plot(gg.iht, gg.iht*0+1, 'k'); hold off;
legend('true h33', 'true h13', 'true h23', 'estim h22', 'estim h13', 'estim h23', 'location', 'southeast');
title('cell 3'); axis tight; set(gca,'ylim',[0,ymax]);

% Print out true and recovered params:
fprintf('\n------------------\n');
fprintf('True stim weights: %.2f, %.2f, %.2f\n', squeeze(gg.k));
fprintf(' Est stim weights: %.2f, %.2f, %.2f\n\n', gg1.k,gg2.k,gg3.k);

% Print out true and recovered params:
fprintf('True dc: %.2f, %.2f, %.2f\n', gg.dc*[1 1 1]);
fprintf(' Est dc: %.2f, %.2f, %.2f\n', gg1.dc,gg2.dc,gg3.dc);
