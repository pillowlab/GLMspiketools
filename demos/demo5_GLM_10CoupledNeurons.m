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

%% Set up basis for coupling filters

% Make basis for self-coupling term
ihbasprs.ncols = 3; % number of basis vectors
ihbasprs.hpeaks = [.001, .005]; % peak of 1st and last basis vector
ihbasprs.b = .001;  % scaling (smaller -> more logarithmic scaling)
ihbasprs.absref = []; % absolute refractory period basis vector (optional)
% Make basis 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dtSp);
nht = length(iht); % number of bins

% Make basis for cross-coupling term
ihbasprs2.ncols = 1;  % number of basis vectors
ihbasprs2.hpeaks = [0.001,.005]; % put peak at 5ms and "effective" 1st peak at 0
ihbasprs2.b = .001;  % smaller -> more logarithmic scaling
ihbasprs2.absref = []; % no abs-refracotry period for this one
% Make basis
[iht2,ihbas2,ihbasis2] = makeBasis_PostSpike(ihbasprs2,dtSp);
nht2 = length(iht2);

% pad to put them into the same time bins
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


%% Set self-coupling filter shapes

wself = [-10; 3.2; -1]; % weights for self-coupling term
ihself = ihbasis*wself; % self-coupling filter
wcpl = 0.5; % weights for cross-coupling term
ihcpl = ihbasis2*wcpl; % cross-coupling filter
clf; plot(iht, exp(ihself), iht, exp(ihcpl), iht, iht*0+1, 'k--');
legend('self-coupling', 'cross-coupling');
xlabel('time lag (s)');
ylabel('gain (sp/s)');


%% Set up multi-neuron GLM

nneur = 10;
k = randn(10,1)*0.5; % stimulus weights
gg.k = permute(k,[2,3,1]);  % stimulus weights
gg.dc = 2+(rand(1,nneur)-0.5);

% Generate random coupling strengths
WW = randn(nneur);  % coupling strengths
WW = WW-diag(diag(WW)-1); % set diagonal to zero (for self-coupling)

gg.iht = iht;
gg.ih = zeros(nht,nneur,nneur);
% Insert coupling filters
for jj = 1:nneur
    gg.ih(:,:,jj) = ihcpl*WW(jj,:); % cross-coupling input weights to neuron jj
    gg.ih(:,jj,jj) = ihself*WW(jj,jj); % self-coupling weights
end

% Plot filters
lw = 2; % linewidth
subplot(131); 
plot(gg.iht, gg.ih(:,:,1),gg.iht,gg.iht*0, 'k--', 'linewidth', lw);
title('incoming filters: cell 1');
subplot(132); 
plot(gg.iht, gg.ih(:,:,2),gg.iht,gg.iht*0, 'k--', 'linewidth', lw);
title('incoming filters: cell 2');
subplot(133); 
plot(gg.iht, gg.ih(:,:,3),gg.iht,gg.iht*0, 'k--', 'linewidth', lw);

%% ===== 2. Run short simulation for visualization purposes ========= %

slen = 200; % Stimulus length (frames) & width (# pixels)
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

slen = 5e5;  % Stimulus length (in bins);  
Stim = stimsd*randn(slen,swid);  % Gaussian white noise stimulus
[tsp,sps,Itot,ispk] = simGLM(gg,Stim);  % run model


%% ===== 4. Do ML fitting for all neurons in a loop =================== %%

% Initialize param struct for fitting 
ggInit = makeFittingStruct_GLM(dtStim,dtSp);  % Initialize params for fitting struct 

% Initialize fields (using h bases computed above)
ggInit.ktbas = 1; % k basis
ggInit.ihbas = ihbas; % h self-coupling basis
ggInit.ihbas2 = ihbas2; % h coupling-filter basis
nktbasis = 1; % number of basis vectors in k basis
nhbasis = size(ihbas,2); % number of basis vectors in h basis
nhbasis2 = size(ihbas2,nneur-1); % number of basis vectors in h basis
ggInit.kt = 1; % initial params from scaled-down sta 
ggInit.k = 1;  % initial setting of k filter
ggInit.ihw = zeros(nhbasis,1); % init params for self-coupling filter
ggInit.ihw2 = zeros(nhbasis2,nneur-1); % init params for cross-coupling filter
ggInit.ih = [ggInit.ihbas*ggInit.ihw ggInit.ihbas2*ggInit.ihw2];
ggInit.iht = iht;
ggInit.dc = 0; % Initialize dc term to zero
ggInit.sps = sps(:,1); % spikes from 1 cell
ggInit.sps2 = sps(:,2:nneur); % spikes from all 
ggInit.couplednums = 2:nneur; % cell numbers of cells coupled to this one 

ggfit(1:nneur) = ggInit; % initialize fitting struct
kfit = zeros(nneur,1);
dcfit = zeros(nneur,1);
for jj = 1:nneur
    couplednums = setdiff(1:nneur,jj);  % cell numbers of cells coupled to this one 

    % Set spike responses for cell 1 and coupled cell
    ggInit.sps = sps(:,jj);
    ggInit.sps2 = sps(:,couplednums);
    ggInit.couplednums = couplednums; % number of cell coupled to this one (for clarity)

    % Do ML fitting
    fprintf('==== Fitting filters to neuron %d ==== \n',  jj)
    opts = {'display', 'iter', 'maxiter', 100};
    ggfit(jj) = MLfit_GLM(ggInit,Stim,opts); % do ML fitting
    kfit(jj) = ggfit(jj).k;
    dcfit(jj) = ggfit(jj).dc;    
end

  
%% ===== 5. Plot true and fitted params for first few cells ============================================= %%

clf reset; colors = get(gca,'colororder');
ncolrs = size(colors,1); % number of traces to plot for each cell 
ymax = max(exp([ggfit(1).ih(:);ggfit(2).ih(:);ggfit(3).ih(:);gg.ih(:)])); % max of y range

for jj = 1:10
    subplot(2,5,jj); % --Spike filters cell 1 % -------------
    ccpl = setdiff(1:nneur,jj); % coupled cells
    cnums = [jj, ccpl];
    plot(gg.iht, exp(gg.ih(:,cnums(1:ncolrs),jj)), ggfit(jj).iht, exp(ggfit(jj).ih(:,1:ncolrs)), '--', 'linewidth', lw);
    hold on; plot(gg.iht, gg.iht*0+1, 'k'); hold off;
    title(sprintf('cell %d filters',jj)); axis tight; set(gca,'ylim',[0,ymax]);
    if jj==1
        ylabel('gain (sp/s)'); xlabel('time after spike (s)');
    end
end

% Print true and recovered params:
fprintf('\n---------------------------------\n');
fprintf('k: true   est   |  dc: true  est\n');
fprintf('---------------------------------\n');
for jj = 1:nneur
    fprintf('  %5.2f  %5.2f        %5.2f %5.2f\n', [gg.k(jj), kfit(jj), gg.dc(jj), dcfit(jj)]);
end
