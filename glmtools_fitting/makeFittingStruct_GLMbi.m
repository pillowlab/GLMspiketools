function gg = makeFittingStruct_GLMbi(krank,varargin)
% gg = makeFittingStruct_GLM(krank,dtStim,dtSp,klength,nkbasis,k0,nhbasis,lasthpeak)
%
% Initialize parameter struct for fitting generalized linear model (GLM),
% with bilinear (i.e., low-rank) parametrization stimulus filter
%
% Inputs:
%       krank = rank of stim filter
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


% ====================================================================
% Set up structure
gg = makeFittingStruct_GLM(varargin{:}); 

% Set additional fields needed by bilinear GLM 
gg.kx = [];
gg.kxbas=[]; 
gg.kbasprs = [];
gg.krank = krank;
gg.ktype = 'bilinear';

% if initial filter passed in, use svd to set up initial K params (bilinearly parametrized)
if (length(varargin)>4) && (~isempty(varargin{3}))
    k0 = varargin{5};
    [u,s,v] = svd(k0);
    mux = mean(k0)';
    nkt = size(gg.ktbas,2);
    nkx = size(k0,2);
    kt = zeros(nkt,krank);
    kx = zeros(nkx,krank);
    for j = 1:krank;
        if v(:,j)'*mux < 0  % Flip sign if shape looks opposite to k0
            u(:,j) = -u(:,j); v(:,j) = -v(:,j); 
        end
        kt(:,j) = (gg.ktbas'*gg.ktbas)\(gg.ktbas'*(u(:,j)*s(j,j)));
        kx(:,j) = v(:,j);
    end
    gg.kxbas = speye(nkx);
    gg.kt = kt;
    gg.kx = kx;
    gg.k = (gg.ktbas*gg.kt)*(gg.kxbas*gg.kx)';
end
