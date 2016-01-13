function gg = makeFittingStruct_GLM_core(sta,DTsim,glmstruct,cellnum)
% gg = makeFittingStruct_GLM_core(sta,DTsim,glmstruct,cellnum)
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
gg.dt = DTsim;
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
ihbasprs.hpeaks = [0 max(DTsim*10,4)];  % Peak location for first and last vectors
ihbasprs.b = 1;  % How nonlinear to make spacings
ihbasprs.absref = []; % absolute refractory period (optional)
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim);
gg.iht = iht;
gg.ihbas = ihbas;
gg.ihbasprs = ihbasprs;
gg.ihw = zeros(size(ihbas,2),1);
gg.ih = gg.ihbas*gg.ihw;  % Set ih to be the current value of ihbas*ihw
        
% % ==================================================================
% If full param struct passed in, match other params as well
if (nargin >= 3) 
    gg.iht = glmstruct.iht;

    % Check for post-spike basis params & update if possible
    if isempty(glmstruct.ih);
        iht = [];
        ih = [];
        ihw = []; ihw2 = [];
        ihbas = []; ihbas2=[];
    else
        if isfield(glmstruct, 'ihbasprs')
            ihbasprs = glmstruct.ihbasprs;
            [iht,ihbas] = makeBasis_PostSpike(ihbasprs,DTsim);
        end
        % Check for separate coupling basis params
        if isfield(glmstruct, 'ihbasprs2');
            ihbasprs2 = glmstruct.ihbasprs2;
            [iht,ihbas2] = makeBasis_PostSpike(ihbasprs2,DTsim,iht);
        else
            ihbas2 = ihbas;
            ihbasprs2 = ihbasprs;
        end
    end
    % Check whether structure is a single-neuron or multi-neuron structure
    ncells = length(glmstruct.dc);
    if (ncells==1)     % Single-neuron
        gg.dc = glmstruct.dc;
        
        gg.ihbas = ihbas;
        gg.ihw = gg.ihbas\glmstruct.ih;
        gg.ihbasprs = ihbasprs;
        gg.ih = gg.ihbas*gg.ihw;  % Set ih to be the current value of ihbas*ihw
        
    else  % multi-neuron structure
        if (nargin == 3)  % Check that cellnum passed in 
            fprintf('proper syntax: gg = makeFittingStruct_GLM_core(sta,DTsim,glmstruct,cellnum)\n');
            error('multi-cell struct passed in, but cell number not specified');
        end

        % Set fields
        gg.dc = glmstruct.dc(cellnum);
        gg.ihbas = ihbas;
        gg.ihbasprs = ihbasprs;
        gg.ihw = gg.ihbas\glmstruct.ih(:,cellnum,cellnum);
        
        gg.couplednums = setdiff(1:ncells,cellnum);
        gg.ihbas2 = ihbas2;
        gg.ihbasprs2 = ihbasprs2;        
        gg.ihw2 = gg.ihbas2\glmstruct.ih(:,gg.couplednums,cellnum);
        gg.ih = [gg.ihbas*gg.ihw, gg.ihbas2*gg.ihw2]; % Set ih to be the current value of ihbas*ihw
    end
end
