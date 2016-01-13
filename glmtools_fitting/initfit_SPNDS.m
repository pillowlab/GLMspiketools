function initfit_SPNDS(gg)
% initfit_SPNDS(gg)
%
% Takes spike times in gg and sets global variables relating to integer
% time indices corresponding to spikes.
%
% Rationale: This frees us from having to pass them between loss functions,
% eg when doing coordinate ascent.
%
% relevant global vars: SPNDS SPNDS2

global SPNDS SPNDS2  % Declare global variables if not passed out

eps = 1e-4; % small number
dt = gg.dt;  % time bin size
 
% Cell's own spike times -----------------------
SPNDS = ceil((gg.tsp-dt*eps)/dt);

% Spike times from coupled cells ---------------
SPNDS2 = [];  
nCoupled = length(gg.tsp2);
for j = 1:nCoupled  
    SPNDS2{j} = ceil((gg.tsp2{j}-dt*eps)/dt);
end

