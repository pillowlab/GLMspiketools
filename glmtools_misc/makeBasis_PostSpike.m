function [iht, ihbas, ihbasis] = makeBasis_PostSpike(ihprs,dt)
% [iht, ihbas, ihbasis] = makeBasis_PostSpike(ihprs,dt)
%
% Make nonlinearly stretched basis consisting of raised cosines
% -------
% Inputs: 
%     prs = param structure with fields:
%            ncols = # of basis vectors
%            hpeaks = 2-vector containg [1st_peak  last_peak], the peak 
%                     location of first and last raised cosine basis vectors
%            b = offset for nonlinear stretching of x axis:  y = log(x+b) 
%                     (larger b -> more nearly linear stretching)
%            absref = absolute refractory period (optional).  If specified,
%                     this param creates an additional "square" basis
%                     vector with support n [0,absref] (i.e., with a hard
%                     cutoff at absref)
%
%     dt = grid of time points for representing basis
%  --------
%  Outputs:  iht = time lattice on which basis is defined
%            ihbas = orthogonalized basis
%            ihbasis = original raised cosine (non-orthogonal) basis 
%
%  -------------
%  Example call:
%
%  ihbasprs.ncols = 5;  
%  ihbasprs.hpeaks = [.1 2];  
%  ihbasprs.b = .5;  
%  ihbasprs.absref = .1;  %% (optional)
%  [iht,ihbas,ihbasis] = makeBasis_PostSpike(ihprs,dt);

ncols = ihprs.ncols;
b = ihprs.b;
hpeaks = ihprs.hpeaks;
if isfield(ihprs, 'absref')
    absref = ihprs.absref;
else
    absref = 0;
end

% Check input values
if (hpeaks(1)+b) < 0 
    error('b + first peak location: must be greater than 0'); 
end
if absref >= dt  % use one fewer "cosine-shaped basis vector
    ncols = ncols-1;
elseif absref > 0
    warning('Refractory period too small for time-bin sizes');
end

% nonlinearity for stretching x axis (and its inverse)
nlin = @(x)log(x+1e-20);
invnl = @(x)exp(x)-1e-20; % inverse nonlinearity
ff = @(x,c,dc)(cos(max(-pi,min(pi,(x-c)*pi/dc/2)))+1)/2; % raised cosine basis vector

% Generate basis of raised cosines
if ncols > 1
    yrnge = nlin(hpeaks+b);        % nonlinearly transformed first & last bumps
    db = diff(yrnge)/(ncols-1);    % spacing between cosine bump peaks
    ctrs = yrnge(1):db:yrnge(2);  % centers (peak locations) for basis vectors
else    
    % if only 1 basis function
    if length(hpeaks)==1
        hpeaks = [0 hpeaks]; % place imaginary 1st bump at 0
    end
    yrnge = nlin(hpeaks+b); % nonlinearly transformed first & last bumps
    ncolseff = 3; % use 3 "effective" number of basis funcs, so left edge decays at 0
    db = diff(yrnge)/(ncolseff-1);    % spacing between cosine bump peaks
    ctrs = yrnge(2);   % centers (peak locations) for basis vectors
end    

% Make basis
mxt = invnl(yrnge(2)+2*db)-b;  % maximum time bin
iht = (dt:dt:mxt)';
nt = length(iht);        % number of points in iht
ihbasis = ff(repmat(nlin(iht+b), 1, ncols), repmat(ctrs, nt, 1), db);

% create first basis vector as step-function for absolute refractory period
if absref >= dt
    ii = find(iht<absref);
    ih0 = zeros(size(ihbasis,1),1);
    ih0(ii) = 1;
    ihbasis(ii,:) = 0;
    ihbasis = [ih0,ihbasis];
end

% compute orthogonalized basis
ihbas = orth(ihbasis);  
