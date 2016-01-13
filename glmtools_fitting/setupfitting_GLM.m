function prs0 = setupfitting_GLM(gg, Stim, maxsize)
%  prs0 = setupfitting_GLM(gg, Stim,maxsize);
%
%  Sets global variables for ML estimation of point process GLM model.
%
%  MSTM = Stimulus filtered by temporal filter basis, 
%  SPNDS = spike indices (vector)
%  SPNDS2 = spike indices of neighor cells (cell array of vectors)
%  MMintrp = sparse matrix for interp of stim to fine time lattice
%  OPTprs = structure with optimization params
%         (ihbas, ihbas2, kxbas, ktbas)
%
%  Inputs: gg = glm param structure
%          Stim = stimulus (time along columns, other dims along rows)
%          maxsize = maximum # floats to store in design matrix 
%                    (Affects # of chunks used to compute likelihood)
%
%  Output: prs0 = initial parameters extracted from gg

% Initialize optimization param structure
global OPTprs
OPTprs = [];  

% ---- Set stimulus filter --------------------------------------------
if strcmp(gg.ktype, 'linear')
    initfit_stimMatrix(gg,Stim); % standard GLM
    ktype = 1;
elseif strcmp(gg.ktype, 'bilinear')
    initfit_stimMatrix_GLMbi(gg,Stim); % bilinearly-parametrized stim filter
    ktype = 2;
elseif strcmp(gg.ktype, 'multiblock');
    ktype = 3;
end

OPTprs.nlfun = gg.nlfun;  % set nonlinearity
initfit_SPNDS(gg);  % set spike indices
initfit_spikeCurrents(gg); % set post-spike current terms
initfit_mask(gg.mask,OPTprs.dt,OPTprs.rlen); % compute mask (time bins to use for likelihood)
initfit_chunks(gg,maxsize); % compute chunks (to avoid out-of-memory errors)

% ---- Extract parameter vector -------------------------------
if (ktype == 1)  % standard GLM
    prs0 = [gg.kt(:); gg.dc; gg.ihw(:); gg.ihw2(:)];

elseif (ktype == 2) % Bilinear k 
    prs0 = [gg.kt(:); gg.kx(:); gg.dc; gg.ihw(:); gg.ihw2(:)];

elseif (ktype == 3)  % mixed rank bilinearly parametrized stim filter
    
%     % Find relevant indices for spatial filters
%     %nkx = OPTprs.nkx;  % fix 
%     xwids = gg.xwids;
%     krank = gg.krank;
%     msk = zeros(nkx,sum(krank));
%     icum = 0;
%     jcum = 0;
%     for jj = 1:length(krank)
%         msk(icum+1:icum+xwids(jj),jcum+1:jcum+krank(jj)) = 1;
%         icum = icum+xwids(jj);
%         jcum = jcum+krank(jj);
%     end
%     ii = find(msk);
% 
%     % Set initial params
%     prs = [gg.kt(:); gg.kx(ii); gg.dc; gg.ihw(:); gg.ihw2(:)];
    
end

