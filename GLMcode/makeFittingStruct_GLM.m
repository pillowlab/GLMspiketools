function gg = makeFittingStruct_GLM(sta,DTsim,glmstruct,cellnumToFit);
% gg = makeFittingStruct_GLM(sta,DTsim,glmstruct,cellnumToFit);
%
% Initialize parameter structure for fitting of GLM model,
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
% set up initial K params
gg.kt = inv(gg.ktbas'*gg.ktbas)*gg.ktbas'*sta;
gg.k = gg.ktbas*gg.kt;

% % ==================================================================
% If full param struct passed in, match other params as well
if (nargin >= 3) 
    gg.dc = glmstruct.dc;
    gg.ih = glmstruct.ih;
    gg.iht = glmstruct.iht;

    %---Extract correct ih basis params, if present----
    if isfield(glmstruct, 'ihbasprs')
        if ~isempty(glmstruct.ihbasprs);
            ihbasprs = glmstruct.ihbasprs;
            [iht,ihbas] = makeBasis_PostSpike(ihbasprs,DTsim);

            % -- Do some error-checking ----
            if length(iht) ~= length(glmstruct.iht)
                error('mismatch between iht and h-kernel params ihbasprs');
            end
            if size(glmstruct.ih,2)>1 & (nargin < 4)
                error('multi-cell glm struct passed in without cell # to fit');
            end

            %--- Put glmstruct params into gg ----
            gg.iht = glmstruct.iht;
            gg.ihbas = ihbas;
            gg.ihbasprs = ihbasprs;
            if nargin == 3  % single-cell only
                gg.ih = inv(ihbas'*ihbas)*ihbas'*glmstruct.ih;
                gg.dc = glmstruct.dc;
            else % mulitcell-cell
                ncells = size(glmstruct.ih,2);
                ih0 = glmstruct.ih(:,:,cellnumToFit);
                ih1 = ih0(:,cellnumToFit);
                ih2 = ih0(:,setdiff(1:ncells,cellnumToFit));
                gg.ih = inv(ihbas'*ihbas)*ihbas'*[ih1 ih2];
                gg.dc = glmstruct.dc(cellnumToFit);
            end

        end
    end


end

