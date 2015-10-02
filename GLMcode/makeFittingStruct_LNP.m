function gg = makeFittingStruct_LNP(sta,DTsim,glmstruct)
% gg = makeFittingStruct_LNP(sta,DTsim,glmstruct,cellnumToFit);
%
% Initialize parameter structure for fitting of LNP model,
% normal parametrization of stim kernel
%
% Inputs:  sta = initial guess at kernel (or use all zeros if unknown)

% Set up structure
gg.k = [];
gg.dc = 0;
gg.ih = [];
gg.iht = [];
gg.ihbas = [];
gg.ihbas2 = [];
gg.ihbasprs = [];
gg.kt = [];
gg.ktbas = [];
gg.kbasprs = [];
gg.tsp = [];
gg.tspi = [];
gg.dt = DTsim;
gg.ih2 = [];
gg.ihbas2 = [];
gg.ihbasprs2 = [];
gg.tsp2 = [];

% === Make temporal basis for stimulus filter =======================
[nkt,nkx] = size(sta);
% % ----- Set up temporal basis for stimulus kernel -----------
kbasprs.neye = 5; % Number of "identity" basis vectors near time of spike;
kbasprs.ncos = 5; % Number of raised-cosine vectors to use  
kbasprs.kpeaks = [0 round(nkt/3)];  % Position of first and last bump (relative to identity bumps)
kbasprs.b = 3; % Offset for nonlinear scaling (larger -> more linear)
ktbas = makeBasis_StimKernel(kbasprs,nkt);
gg.ktbas = ktbas;
gg.kbasprs = kbasprs;

% % ==================================================================
% set up initial K params
gg.kt = inv(gg.ktbas'*gg.ktbas)*gg.ktbas'*sta;
gg.k = gg.ktbas*gg.kt;

% % ==================================================================
% If full param struct passed in, match other params as well
if (nargin >= 3) 
    gg.dc = glmstruct.dc;
end

