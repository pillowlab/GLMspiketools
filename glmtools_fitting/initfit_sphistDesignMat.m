function Xstruct = initfit_sphistDesignMat(gg,Xstruct)
%  Xstruct = initfit_sphistDesignMat(gg,Xstruct)
%
%  Sets parameters relating to optimization of
%  spike-filter terms in the point process GLM and inserts them into design
%  matrix structure 'Xstruct'


% Determine number of parameters from self-coupling filter
if ~isempty(gg.ihw), 
    nh=size(gg.ihbas,2);
else
    nh=0;   % Number of basis vectors describing ih
end

% Determine number of coupled neurons and parameters from coupling filters
nCoupled = length(gg.couplednums); % Number of coupled cells
if ~isempty(gg.ihw2)  
    nh2 = size(gg.ihbas2,2);  % # params describing coupling filters
else  % No coupling terms
    nh2=0;
end

% Set flag for whether or not to worry about any of ih params
if (nh2 == 0) && (nh == 0)
    Xstruct.ihflag = 0;
else
    Xstruct.ihflag = 1;
end

% Compute binned spike times 
dt = gg.dtSp;  % time bin size
Xstruct.spInds = find(gg.sps);

% Create space for design matrix
Xsphist = zeros(Xstruct.rlen,nh+nh2*nCoupled);

% Create design matrix columns for neuron's own spike-history
twin = [1 Xstruct.rlen]; % time window for convolution (entire length)
if nh>0
    Xsphist(:,1:nh) = spikefilt_mex(Xstruct.spInds, gg.ihbas, twin);
end

% Create design matrix columns for input from coupled neurons
for jcpl = 1:nCoupled
    spInds_jcpl = find(gg.sps2(:,jcpl)); % spike times of coupled neuron
    inds = nh+nh2*(jcpl-1)+1:nh+nh2*jcpl; % column indices
    Xsphist(:,inds) = spikefilt_mex(spInds_jcpl,gg.ibas2,twin);
end

% Add to Xstruct
Xstruct.Xsphist = Xsphist;