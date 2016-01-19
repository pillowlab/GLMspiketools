function [prs0,Xstruct] = setupfitting_GLM(gg, Stim)
%  prs0 = setupfitting_GLM(gg, Stim);
%
%  Sets global variables for ML estimation of point process GLM model.
%
%  MSTM = Stimulus filtered by temporal filter basis, 
%  SPNDS = spike indices (vector)
%  SPNDS2 = spike indices of neighor cells (cell array of vectors)
%  OPTprs = structure with optimization params
%         (ihbas, ihbas2, kxbas, ktbas)
%
%  Inputs: 
%     gg = glm param structure
%     Stim = stimulus (time along columns, other dims along rows)
%     maxsize = maximum # floats to store in design matrix 
%                     (Affects # of chunks used to compute likelihood)
%  Output: 
%     prs0 = initial parameters extracted from gg
%     Regressors = struct with design matrix for stim and spike-history terms

%
% Initialize optimization param structure

% ---- Create struct and make stimulus design matrix ---------------------
if strcmp(gg.ktype, 'linear') % standard GLM
    Xstruct = initfit_stimDesignMat(gg,Stim);  % create design matrix structure
    
    % extract parameter vector
    prs0 = [gg.kt(:); gg.dc; gg.ihw(:); gg.ihw2(:)];

elseif strcmp(gg.ktype, 'bilinear') % bilinearly-parametrized stim filter
    Xstruct = initfit_stimMatrix_GLMbi(gg,Stim); % create design matrix structure

    % extract parameter vector    
    prs0 = [gg.kt(:); gg.kx(:); gg.dc; gg.ihw(:); gg.ihw2(:)];

else
    error('unknown filter type (allowed types are ''linear'' or ''bilinear'')');
end

% ---- Make spike-history design matrix -----------------------------------
Xstruct = initfit_sphistDesignMat(gg,Xstruct); 



% set nonlinearity
Xstruct.nlfun = gg.nlfun;  

% set spike indices
eps = 1e-6; % small number to make sure spikes not on a bin edge
dt = gg.dtSp;  % time bin size
 
% Cell's own spike times -----------------------
spInds = ceil((gg.tsp-dt*eps)/dt);

% Spike times from coupled cells ---------------
nCoupled = length(gg.tsp2);
spInds2 = cell(nCoupled,1);
for j = 1:nCoupled  
    spInds2{j} = ceil((gg.tsp2{j}-dt*eps)/dt);
end


% compute mask (time bins to use for likelihood)
initfit_mask(gg.mask,Xstruct.dt,Xstruct.rlen); 


