function setSPNDS(gg);
% setSPNDS(gg);
%
% 1. Map scalar spike times into integer sampling lattice
% 2. Compute matrix for mapping h into correct sampling
%
% Sets global variables: SPNDS, SPNDS2, MSP

global SPNDS SPNDS2 MSP

% ---  Compute spike times as integers in new time sampling lattice -----
dt = gg.dt;
gg.tsp = gg.tsp(:);
SPNDS = ceil((gg.tsp-dt*.001)/dt); % own spike train
SPNDS2 = [];  % coupled spike trains
tsp2 = gg.tsp2;
nCoupled = length(tsp2);
for j = 1:nCoupled  % spike trains of other cells
    SPNDS2{j} = ceil((tsp2{j}-dt*.001)/dt);
end

% Compute MSP matrix
if ~isempty(gg.iht)
    ihdt = diff(gg.iht(1:2));
    hlen = length(gg.iht);
    if (round(ihdt/dt) == ihdt/dt)
        MSP = makeInterpMatrix(hlen,ihdt/dt);
        MSP(1:ihdt/dt-1,:) = [];  % Remove bins that lie prior to time 0
        MSP(1,1) = 0;  % set first entry to 0, so as not to inject during spike time
    else
        error(['setSPNDS: time binning iht must be multiple of gg.dt\n']);
    end
end
