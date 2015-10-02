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
%        'fnlin' - nonlinearity (exponential by default)
%        'dc' - dc input to cell
%        'ih' - post-spike current
%        'ihbas' - basis for post-spike current
%        'iht' - time lattice for post-spike current
%        'dtsim' - default time bin size for simulation
%        'kbasprs' - basis for stim filter 
%        'ihbasprs' - basis for post-spike current

% Create a default (temporal) stimulus filter
tk = [0:nkt-1]';
b1 = nkt/32; b2 = nkt/16;
k1 = 1/(gamma(6)*b1)*(tk/b1).^5 .* exp(-tk/b1);  % Gamma pdfn
k2 = 1/(gamma(6)*b2)*(tk/b2).^5 .* exp(-tk/b2);  % Gamma pdf
k = flipud(k1-k2./1.5);

% -- Nonlinearity -------
nlinF = @exp;
% Other natural choice:  nlinF = @(x)log(1+exp(x));

% --- Make basis for post-spike (h) current ------
ihbasprs.ncols = 5;  % Number of basis vectors for post-spike kernel
ihbasprs.hpeaks = [.1 2];  % Peak location for first and last vectors
ihbasprs.b = .5;  % How nonlinear to make spacings
ihbasprs.absref = .1; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,dt);
ih = ihbasis*[-10 -5 0 2 -2]';  % h current

% Place parameters in structure
S = struct(...
    'type', 'glm', ...
    'k', 2*k, ... % stimulus filter
    'nlfun', nlinF, ...  % nonlinearity
    'dc', 3, ...         % dc current 
    'ih', ih, ...        % post-spike current
    'iht', iht, ...      % time indices of aftercurrent
    'dt', dt, ...
    'kbasprs', [], ... % params for stim filter basis
    'ihbasprs', ihbasprs);  % params for ih basis
