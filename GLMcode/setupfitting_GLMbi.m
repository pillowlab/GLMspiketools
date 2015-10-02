function setupfitting_GLMbi(gg, Stim, maxsize);
%  setupfitting_GLMbi(gg, Stim,maxsize);
%
%  Sets global variables for ML estimation of GLM model using 'maxli_GLM'
%
%  MSTM = Stimulus filtered by temporal filter basis, 
%  MSP = Interpolates ih current to fine temporal sampling
%  SPNDS = spike indices (vector)
%  SPNDS2 = spike indices of neighor cells (cell array of vectors)
%  MMintrp = sparse matrix for interp of stim to fine time lattice
%  OPRS = structure with optimization params
%         (ihbas, ihbas2, kxbas, ktbas)

global SPNDS SPNDS2 MSP MSTM  MMntrp OPRS

setSPNDS(gg);  % set SPNDS and SPNDS2 (integer spike times in time lattice)

% ---- Set post-spike filter params -------
OPRS = [];
OPRS.dt = gg.dt;
OPRS.iht = gg.iht;
if isempty(gg.ih), 
     OPRS.ihflag=0; 
     OPRS.nh=0;
else
    OPRS.ihflag=1;  
    OPRS.nh=size(gg.ihbas,2);
    OPRS.ihbas= gg.ihbas;

end
nCoupled = length(SPNDS2);
if isempty(gg.ihbas2)
    OPRS.ih2flag=0; 
    OPRS.nh2=0;
else
    OPRS.ih2flag = 1; 
    OPRS.nh2 = size(gg.ihbas2,2);
    OPRS.ihbas2 = gg.ihbas2;
end

% ---- Set up stim processing params ----------
krank = size(gg.kt,2);
OPRS.nx = size(gg.kxbas,2);
OPRS.nkx = krank*OPRS.nx;
OPRS.nt = size(gg.ktbas,2);
OPRS.nkt = krank*OPRS.nt;
OPRS.krank = krank;
OPRS.ktbas = gg.ktbas;
OPRS.kxbas = gg.kxbas;

% ---- Filter stimulus with spatial and temporal bases -----
[slen,swid] = size(Stim);
rlen = round(slen/OPRS.dt);
xfltStim = Stim*gg.kxbas;
nx = OPRS.nx; nt = OPRS.nt;
ncols = nx*nt;
MSTM = zeros(slen,ncols);
for i = 1:nx
    for j = 1:nt
        MSTM(:,(i-1)*nt+j) = sameconv(xfltStim(:,i),gg.ktbas(:,j));
    end
end

% ------ Compute Spike matrix MSP - maps ih into correct time sampling ---
dt = OPRS.dt;
if OPRS.ihflag
    OPRS.ihbas = MSP*OPRS.ihbas; % filter ih basis
end
if OPRS.ih2flag
    OPRS.ihbas2 = MSP*OPRS.ihbas2; % filter ih2 basis
end
OPRS.tspi = gg.tspi(1);
OPRS.ncpl = nCoupled;

% --- Set up chunking in time, for managing memory  --------
if isempty(OPRS.tspi),
    istrt1 = 0;
elseif (OPRS.tspi(1) == 1),
    istrt1 = 0;
else,                   
    istrt1 = SPNDS(OPRS.tspi(1)-1);
end
istrt2 = ceil(istrt1*dt)/dt;
istrt = unique([istrt1, istrt2]);  % first two indices

rlenLi = rlen-istrt2;  % Total length to be chunked up -------
slenLi = slen-istrt2*dt;

ndt = round(1/OPRS.dt);
nhtot = OPRS.nh+OPRS.nh2*nCoupled;
totNums = (nhtot+1)*rlenLi + size(MSTM,2)*slenLi;
nchunks = ceil(totNums/maxsize);
dchk0 = floor(slenLi/nchunks);
dchunk = dchk0*ndt; 

% ichunk are chunks in fine time; jchunk are coarse time
ichunk = [istrt, (istrt2+dchunk):dchunk:rlen];
if ichunk(end)~=rlen
    ichunk(end+1) = rlen;
end
nchunks = length(ichunk)-1;
jchunk = [floor(ichunk(1:end-1)'*dt)-1, ceil(ichunk(2:end)'*dt)+1];
jchunk(1,1) = max(jchunk(1,1),0);
jchunk(nchunks,2) = min(slen,jchunk(nchunks,2));

OPRS.nchunks = nchunks;
OPRS.ichunk = ichunk;
OPRS.jchunk = jchunk;

% ----- Make time interpolation matrix for stimulus -----------
MMntrp = makeInterpMatrix2(max(diff(jchunk')),ndt); 

% ----- Nonlinearity -----------
if isfield(gg, 'nlfun')
    OPRS.nlfun = gg.nlfun;
else
    OPRS.nlfun = @expfun;
end
