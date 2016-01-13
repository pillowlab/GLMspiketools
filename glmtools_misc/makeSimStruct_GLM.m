function S = makeSimStruct_GLM(nkt,dt)
%  S = makeSimStruct_GLM(nkt,dt);
%
%  Creates a structure with default parameters for a GLM model
%
%  Input: nkt = number of time bins for stimulus filter
%         dt = time binning for sampling post-spike kernel ('ih')
%
%  Struct fields (model params):  
%        'filt' - stimulus filter
%        'nlfun' - nonlinearity (exponential by default)
%        'dc' - dc input to cell
%        'ih' - post-spike current
%        'ihbas' - basis for post-spike current
%        'iht' - time lattice for post-spike current
%        'dtsim' - default time bin size for simulation
%        'ihbasprs' - basis for post-spike current

% Create a default (temporal) stimulus filter
tk = (0:nkt-1)'; % time indices for made-up filter 
b1 = nkt/32; b2 = nkt/16;
k1 = 1/(gamma(6)*b1)*(tk/b1).^5 .* exp(-tk/b1);  % Gamma pdfn
k2 = 1/(gamma(6)*b2)*(tk/b2).^5 .* exp(-tk/b2);  % Gamma pdf
k = flipud(k1-k2./1.5);
k = k./norm(k)/2;

% % ----- Represent this filter (approximately) in temporal basis -------
ktbasprs.neye = min(5,nkt); % Number of "identity" basis vectors near time of spike;
ktbasprs.ncos = min(5,nkt); % Number of raised-cosine vectors to use  
ktbasprs.kpeaks = [0 ((nkt-ktbasprs.neye)/2)];  % Position of 1st and last bump 
ktbasprs.b = 1; % Offset for nonlinear scaling (larger -> more linear)
ktbas = makeBasis_StimKernel(ktbasprs,nkt);
k = ktbas*(ktbas\k);

% -- Nonlinearity ------- 
nlinF = @expfun;
% Other natural choice:  nlinF = @(x)log(1+exp(x));

% --- Make basis for post-spike (h) filter ------
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [0 max(dt*10,4)];  % Peak location for first and last vectors
ihbasprs.b = 1;  % How nonlinear to make spacings
ihbasprs.absref = []; % absolute refractory period (optional)
[iht,~,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
ih = ihbasis*[-10 .5 1 -.5 -.5]';  % h current

% Place parameters in structure
S = struct(...
    'k', 2*k, ... % stimulus filter
    'nlfun', nlinF, ...  % nonlinearity
    'dc', 3, ...         % dc input (constant) 
    'ih', ih, ...        % post-spike current
    'iht', iht, ...      % time indices of aftercurrent
    'dt', dt, ...
    'ihbasprs', ihbasprs, ...  % params for ih basis
    'ktbasprs', ktbasprs);     % params for k time-basis