function initfit_spikeCurrents(gg)
%  initfit_spikeCurrents(gg)
%
%  Sets global variables relating to optimization of
%  spike-filter terms in the point process GLM
%
%  OPTprs = global variable (structure) for optimization params


global OPTprs

OPTprs.iht = gg.iht;

% --- Recurrent (self) post-spike filter -----------------
if ~isempty(gg.ihw), 
    OPTprs.nh=size(gg.ihbas,2);
    OPTprs.ihbas= interp_spikeFilt(gg.ihbas, gg.iht, gg.dt);
else
    OPTprs.nh=0;  % Number of basis vectors describing ih
end

% --- Coupling Filters -----------------------------------
OPTprs.ncpl = length(gg.couplednums); % Number of coupled cells

if ~isempty(gg.ihw2)  
    OPTprs.nh2 = size(gg.ihbas2,2);  % # params describing coupling filters
    OPTprs.ihbas2 = interp_spikeFilt(gg.ihbas2, gg.iht, gg.dt); % relevant basis
else  % No coupling terms
    OPTprs.nh2=0;
end

% Flag for whether or not to worry about ih params
if (OPTprs.nh2 == 0) && (OPTprs.nh == 0)
    OPTprs.ihflag = 0;
else
    OPTprs.ihflag = 1;
end
