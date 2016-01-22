function gg = makeFittingStruct_GLMbi(sta,DTsim,krank,varargin)
% gg = makeFittingStruct_GLMbi(sta,DTsim,krank,glmstruct);
%
% Initialize parameter structure for fitting of GLM model
% with bilinear parametrization of stim kernel
%
% Inputs:  sta = initial guess at kernel (or use all zeros if unknown)
%          DTsim = bin size for simulations / likelihood calculation
%          krank = rank of stim kernel "k"
%          glmstruct (optional) = glm sim structure to extract from


% ====================================================================
% Set up structure
gg = makeFittingStruct_GLM(sta,DTsim,varargin{:}); 
gg.kx = [];
gg.kxbas=[]; 
gg.kbasprs = [];
gg.krank = krank;

% % ==================================================================
% use svd to set up initial K params (bilinearly parametrized)
[u,s,v] = svd(sta);
mux = mean(sta)';
nkt = size(gg.ktbas,2);
nkx = size(sta,2);

kt = zeros(nkt,krank);
kx = zeros(nkx,krank);
for j = 1:krank;
    if v(:,j)'*mux < 0  % Flip sign if match is better
        u(:,j) = -u(:,j); v(:,j) = -v(:,j);
    end
    kt(:,j) = (gg.ktbas'*gg.ktbas)\(gg.ktbas'*(u(:,j)*s(j,j)));
    kx(:,j) = v(:,j);
end
gg.kxbas = speye(nkx);
gg.kt = kt;
gg.kx = kx;
gg.k = (gg.ktbas*gg.kt)*(gg.kxbas*gg.kx)';

gg.ktype = 'bilinear';