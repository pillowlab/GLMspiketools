function gg = makeFittingStruct_GLM(sta,DTsim,varargin)
% gg = makeFittingStruct_GLM(sta,DTsim,glmstruct,cellnumToFit);
%
% Initialize parameter structure for fitting of GLM model,
% normal parametrization of stim kernel
%
% Inputs:  sta = initial guess at kernel (or use all zeros if unknown)
%          DTsim = bin size for simulations / likelihood calculation
%          glmstruct (optional) = glm sim structure to extract from


% ====================================================================
% Set up structure
gg = makeFittingStruct_GLM_core(sta,DTsim,varargin{:}); 

% % ==================================================================
% set up initial K params (stim filter, linearly parametrized)
gg.kt = inv(gg.ktbas'*gg.ktbas)*gg.ktbas'*sta;
gg.k = gg.ktbas*gg.kt;

gg.ktype = 'linear';

