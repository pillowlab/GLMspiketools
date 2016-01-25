function Xstruct = initfit_sphistDesignMat(gg,Xstruct)
%  Xstruct = initfit_sphistDesignMat(gg,Xstruct)
%
%  Sets parameters relating to optimization of
%  spike-filter terms in the point process GLM and inserts them into design
%  matrix structure 'Xstruct'


% Determine # of parameters for self-coupling filter
if ~isempty(gg.ihw),  nh=size(gg.ihbas,2);
else nh=0;  
end

% Determine # of params for coupling filters
if ~isempty(gg.ihw2), nh2 = size(gg.ihbas2,2); 
else  nh2=0;
end

% Determine # of coupled neurons
nCoupled = length(gg.couplednums); % Number of coupled cells

% Set flag for existence of coupling filters
if (nh+nh2 == 0),  Xstruct.ihflag = false;
else Xstruct.ihflag = 1;
end

% Vector of binary (0 or 1) spike counts
Xstruct.bsps = sparse(gg.sps>0);

% ---- Create Spike-history Design Matrix ------------------------
Xsp = zeros(Xstruct.rlen,nh+nh2*nCoupled); % allocate

% Design matrix for self-coupling filter 
if nh>0
    Xsp(:,1:nh) = spikefilt(full(double(gg.sps>0)),gg.ihbas);
end

% Design matrix for cross-coupling filters
for jcpl = 1:nCoupled
    inds = nh+nh2*(jcpl-1)+1:nh+nh2*jcpl; % column indices
    Xsp(:,inds) = spikefilt(full(double(gg.sps2(:,jcpl))),gg.ihbas2); % insert into design matrix
end

% ---- Set fields of Xstruct -------------------------------------
Xstruct.nh = nh;
Xstruct.nh2 = nh2;
Xstruct.nCoupled = nCoupled;
Xstruct.Xsp = Xsp;