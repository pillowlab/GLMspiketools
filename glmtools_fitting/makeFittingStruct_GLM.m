function gg = makeFittingStruct_GLM(dtStim,dtSp,klength,nkbasis,k0,nhbasis,lasthpeak)
% gg = makeFittingStruct_GLM(dtStim,dtSp,klength,nkbasis,k0,nhbasis,lasthpeak)
%
% Initialize parameter structure for fitting GLM,
% with normal parametrization of stim kernel
%
% Inputs:
%      dtStim = bin size of stimulus (in s)
%        dtSp = bin size for spike train (in s)
%     klength = temporal length of stimulus filter, in # of bins (optional)
%     nkbasis = # temporal basis vectors for stim filter         (optional)
%     nhbasis = # temporal basis vectors for spike hist filter h (optional)
%   lasthpeak = time of peak of last h basis vector, in s        (optional)
%          k0 = initial estimate of stimulus filter              (optional)   
%
% Outputs:
%           gg = GLM-fitting param structure

% ---- Set up fitting structure -------------------------------
gg.k = [];  % Actual stim filter k
gg.ih = []; % Actual post-spike filter ih
gg.dc = 0;  % "dc" or constant input (determines baseline spike rate)
gg.ihw = [];  % h filter weights 
gg.iht = [];  % h filter time points
gg.ihbas = []; % basis for h filter
gg.ihbasprs = []; % parameters for basis for h-filter
gg.kt = [];       % basis weights for stimulus filter k
gg.ktbas = [];    % temporal basis for stimulus filter k
gg.ktbasprs = []; % parameters for basis for k-filter
gg.nlfun = @expfun; % default nonlinearity: exponential
gg.sps = [];  % spike times (in s)
gg.mask = []; % list of intervals to ignore when computing likelihood
gg.dtStim = dtStim;  % time bin size for stimulus 
gg.dtSp = dtSp;      % time bin for spike train 
gg.ihw2 = [];      % weights for coupling filters 
gg.ihbas2 = [];    % basis for coupling filters
gg.ihbasprs2 = []; % parameters for coupling filter basis 
gg.sps2 = [];      % spike times of coupled neurons
gg.couplednums = []; % numbers of coupled cells

% % ----- Set up temporal basis for stimulus kernel -----------
if nargin > 2
    assert((klength>nkbasis), 'klength should be bigger than number of temporal basis vectors');

    ktbasprs.neye = 0; % number of "identity" basis vectors
    ktbasprs.ncos = nkbasis; % Number of raised-cosine vectors to use
    ktbasprs.kpeaks = [0 klength*(1 - 1.5/nkbasis)];  % Position of 1st and last bump
    ktbasprs.b = 10; % Offset for nonlinear scaling (larger -> more linear)
    [~,ktbasis] = makeBasis_StimKernel(ktbasprs,klength);
    gg.ktbas = ktbasis;
    gg.ktbasprs = ktbasprs;
    
end

if (nargin > 4) && (~isempty(k0))
    % initialize k filter in this basis
    gg.kt = (gg.ktbas'*gg.ktbas)\(gg.ktbas'*k0);
    gg.k = gg.ktbas*gg.kt;
end

% ----- Set up basis for post-spike filter -----------------------
if nargin > 5
    ihbasprs.ncols = nhbasis;  % Number of basis vectors for post-spike kernel
    ihbasprs.hpeaks = [dtSp lasthpeak];  % Peak location for first and last vectors
    ihbasprs.b = lasthpeak/5;  % How nonlinear to make spacings (rough heuristic)
    ihbasprs.absref = []; % absolute refractory period (optional)
    [iht,ihbas] = makeBasis_PostSpike(ihbasprs,dtSp);
    gg.iht = iht;
    gg.ihbas = ihbas;
    gg.ihbasprs = ihbasprs;
    gg.ihw = zeros(size(ihbas,2),1);
    gg.ih = gg.ihbas*gg.ihw;  % Set ih to be the current value of ihbas*ihw
end

% specify type of parametrization of filter ('linear' vs. 'bilinear')
gg.ktype = 'linear';

