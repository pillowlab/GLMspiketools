% Unit test to check the consistency of conditional intensity computed
% during simulation and during fitting, using spikes from generated from a
% simulation of coupled neurons in demo3  


% 1. First, run relevant section of demo3

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


% Run short simulation for visualization purposes

slen = 50; % Stimulus length (frames) & width (# pixels)
swid = size(ggsim.k,2); % stimulus width
Stim = 2*randn(slen,swid);  % Gaussian white noise stimulus
[tsp,sps,Itot,Istm] = simGLM(ggsim, Stim); % Simulate GLM response
Isp = Itot-Istm; % net spike-history output


%% 2. compute conditional intensity again during log-likelihood evaluation

%  Make param object with "true" params
ggTrue = makeFittingStruct_GLM(dtStim,dtSp);
ggTrue.k = ggsim.k(:,:,1); % % insert stimulus filter
ggTrue.dc = ggsim.dc(1); % insert dc value
ggTrue.ih = ggsim.ih(:,:,1); % insert spike-history filter

ggTrue.sps = sps(:,1);  % spike times
ggTrue.sps2 = sps(:,2); % spike times of coupled neuron
ggTrue.couplednums = 2; % number of cell coupled to this one (for clarity)

[logliTrue2,rrT2,tt2,Itot2,Istm2,Ih2] = neglogli_GLM(ggTrue,Stim);

% ---- Unit tests -------

% On total log-conditional intensity
assert(max(abs(Itot(:,1)-Itot2))<1e-8,'Unit failed: log-conditional intensity consistency');

% On just the stimulus component of the conditional intensity
assert(max(abs(Istm(:,1)-Istm2))<1e-8,'Unit failed: spike-history filter output consistency');

%% 3. Do same for other neuron

%  Make param object with "true" params
ggTrue = makeFittingStruct_GLM(dtStim,dtSp);
ggTrue.k = ggsim.k(:,:,2); % % insert stimulus filter
ggTrue.dc = ggsim.dc(2); % insert dc value
ggTrue.ih = fliplr(ggsim.ih(:,:,2)); % insert spike-history filter

ggTrue.sps = sps(:,2);  % spike times
ggTrue.sps2 = sps(:,1); % spike times of coupled neuron
ggTrue.couplednums = 1; % number of cell coupled to this one (for clarity)

[logliTrue3,rrT3,tt3,Itot3,Istm3,Ih3] = neglogli_GLM(ggTrue,Stim);

% ---- Unit tests -------

% On total log-conditional intensity
assert(max(abs(Itot(:,2)-Itot3))<1e-8,'Unit failed: log-conditional intensity consistency');

% On just the stimulus component of the conditional intensity
assert(max(abs(Istm(:,2)-Istm3))<1e-8,'Unit failed: spike-history filter output consistency');

% =================================== 
% Report if passed
% ===================================
fprintf('Unit passed: condIntensityConsistency_cpl\n');

