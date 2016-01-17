function gg = makeFittingStruct_GLM_core(sta,dtStim,dtSp)
% gg = makeFittingStruct_GLM_core(sta,dtStim,dtSp)
%
%  Core common to both makeFittingStruct_GLM and makeFittingStruct_GLMbi
%
% Inputs:  sta = initial guess at kernel (or use all zeros if unknown)
%          DTsim = bin size for simulations / likelihood calculation
%          glmstruct (optional) = glm sim structure to extract from


% Set up structure =====================================
gg.k = [];  % Actual stim filter k
gg.ih = []; % Actual post-spike filter ih
gg.dc = 0;
gg.ihw = [];  % ih weights 
gg.iht = [];  % ih time points
gg.ihbas = [];
gg.ihbas2 = [];
gg.ihbasprs = [];
gg.kt = [];
gg.ktbas = [];
gg.ktbasprs = [];
gg.nlfun = @expfun; % default nonlinearity: exponential
gg.tsp = [];
gg.mask = [];
gg.dtStim = dtStim;
gg.dtSp = dtSp;
gg.ihw2 = [];
gg.ihbas2 = [];
gg.ihbasprs2 = [];
gg.tsp2 = [];
gg.couplednums = [];

nkt = size(sta,1);

% % ----- Set up temporal basis for stimulus kernel -----------
ktbasprs.neye = min(5,nkt); % Number of "identity" basis vectors near time of spike;
ktbasprs.ncos = min(5,nkt); % Number of raised-cosine vectors to use  
ktbasprs.kpeaks = [0 ((nkt-ktbasprs.neye)/2)];  % Position of 1st and last bump 
ktbasprs.b = 1; % Offset for nonlinear scaling (larger -> more linear)
ktbas = makeBasis_StimKernel(ktbasprs,nkt);
gg.ktbas = ktbas;
gg.ktbasprs = ktbasprs;

% ======================================================================
% Make default basis for post-spike filter

ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [0 max(dtSp*10,4)];  % Peak location for first and last vectors
ihbasprs.b = 1;  % How nonlinear to make spacings
ihbasprs.absref = []; % absolute refractory period (optional)
[iht,ihbas] = makeBasis_PostSpike(ihbasprs,dtSp);
gg.iht = iht;
gg.ihbas = ihbas;
gg.ihbasprs = ihbasprs;
gg.ihw = zeros(size(ihbas,2),1);
gg.ih = gg.ihbas*gg.ihw;  % Set ih to be the current value of ihbas*ihw
        
