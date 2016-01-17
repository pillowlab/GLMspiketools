function [spInds,spInds2] = initfit_spikeInds(gg)
% [spInds,spInds2] = initfit_spikeInds(gg)
%
% Takes spike times in gg.tsp and gg.tsp_cpl and compute indices
% in discrete time bins of width gg.dtSp.


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

