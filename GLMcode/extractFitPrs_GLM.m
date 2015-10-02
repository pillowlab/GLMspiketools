function prs = extractFitPrs_GLM(gg,Stim,MAXSIZE);
% prs = extractFitPrs_GLM(gg,Stim,MAXSIZE);
%
% Set global variables for fitting and extract the 
% the parameters needed for fitting the (vanilla) GLM model

setupfitting_GLM(gg,Stim,MAXSIZE);  % Precompute quantities for optimization

prs = [gg.kt(:); gg.dc; gg.ih(:); gg.ih2(:)];

