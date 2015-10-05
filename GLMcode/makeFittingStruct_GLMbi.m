function gg = makeFittingStruct_GLMbi(sta,krank,DTsim,glmstruct);
% gg = makeFittingStruct_GLMbi(sta,krank,DTsim,glmstruct);
%
% Initialize parameter structure for fitting of GLM model
% with bilinear parametrization of stim kernel
%
% Inputs:  krank = rank of stim kernel "k"
%          sta = initial guess at kernel (or use all zeros if unknown)

% Set up structure
gg.k = [];
gg.dc = 0;
gg.ih = [];
gg.iht = [];
gg.ihbas = [];
gg.ihbas2 = [];
gg.ihbasprs = [];
gg.tsp = [];
gg.tspi = [];
gg.dt = DTsim;
gg.ktbas=[];
gg.kxbas=[]; 
gg.kbasprs = [];
gg.kt = [];
gg.kx = [];
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
gg.kbasprs = kbasprs;

% ======================================================================
% Set up basis for post-spike kernel

ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [DTsim*10 2];  % Peak location for first and last vectors
ihbasprs.b = .4;  % How nonlinear to make spacings
ihbasprs.absref = DTsim*10; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim);
gg.iht = iht;
gg.ihbas = ihbas;
gg.ihbasprs = ihbasprs;
gg.ih = zeros(size(ihbas,2),1);

% % ==================================================================
% use svd to set up initial K params
[u,s,v] = svd(sta);
mux = mean(sta)';
nkt = size(ktbas,2);

kt = zeros(nkt,krank);
kx = zeros(nkx,krank);
for j = 1:krank;
    if v(:,j)'*mux < 0  % Flip sign if match is better
        u(:,j) = -u(:,j); v(:,j) = -v(:,j);
    end
    kt(:,j) = (ktbas'*ktbas)\(ktbas'*(u(:,j)*s(j,j)));
    kx(:,j) = v(:,j);
end
gg.ktbas = ktbas;
gg.kxbas = eye(nkx);
gg.kt = kt;
gg.kx = kx;
gg.k = (gg.ktbas*gg.kt)*(gg.kxbas*gg.kx)';

% % ==================================================================
% If full param struct passed in, match post-spike kernel as well
if nargin > 3
    ihbasprs = glmstruct.ihbasprs;
    [iht,ihbas] = makeBasis_PostSpike(ihbasprs,DTsim);
    if isfield(glmstruct, 'ihbas')
        ih = glmstruct.ihbas*glmstruct.ih;
    else
        ih = glmstruct.ih;
    end
    gg.iht = iht;
    gg.ihbas = ihbas;
    gg.ih = (ihbas'*ihbas)\(ihbas'*ih);
    gg.dc = glmstruct.dc;
end


