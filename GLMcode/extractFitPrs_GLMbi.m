function prs = extractFitPrs_GLMbi(gg,Stim,MAXSIZE);
% prs = extractFitPrs_GLMbi(gg,Stim,MAXSIZE););
%
% Set global variables for fitting and extract the 
% the parameters needed for fitting the (vanilla) GLM model


setupfitting_GLMbi(gg,Stim,MAXSIZE);  % Precompute quantities for optimization

% Set initial params
prs = [gg.kt(:); gg.kx(:); gg.dc; gg.ih(:); gg.ih2(:)];

