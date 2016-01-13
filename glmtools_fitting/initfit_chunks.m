function initfit_chunks(gg,maxsize)
% initfit_chunks(gg,maxsize)
%
% Set up the chunks for computing likelihood without running out of memory
%
% Note: if mask set in OPTprs has multiple windows, uses the mask to define
% chunks (i.e., won't break long blocks to avoid memory errors unless these
% blocks appear in mask).

global OPTprs MSTM

dt = gg.dt;  % Time bin size
ndt = round(1/dt); % # fine time bins per coarse timebin
lastInd = OPTprs.iiLi(end);


if size(OPTprs.mask,1) > 1  % Use mask itself for chunking

    nchunks = size(OPTprs.mask,1);
    ichunk = OPTprs.mask;
    jchunk = round(OPTprs.mask*dt);
    jchunk(:,1) = max(1,jchunk(:,1)-1);
    jchunk(:,2) = min(OPTprs.slen,jchunk(:,2)+1);

else  % Chunk to avoid out-of-memory errors

    ilims = [OPTprs.iiLi(1)-1, lastInd]; % fine index range (from mask)
    rlenLi = diff(ilims);  % total length to be chunked up  (fine timescale)
    slenLi = ceil(rlenLi*dt);  % total length to be chunked up  (coarse timescale)

    nhprs = OPTprs.nh+OPTprs.nh2*OPTprs.ncpl; % # params for h currents
    ntotprs = nhprs*rlenLi + size(MSTM,2)*slenLi; % total # params

    % Compute number of chunks
    nchunks = ceil(ntotprs/maxsize);
    dchk0 = floor(slenLi/nchunks);
    dchunk = dchk0*ndt;

    % ====  ichunk: chunks in fine time =====================================
    i2 = ceil(ilims(1)*dt)/dt;  % start of second chunk (break in coarse grid)
    initialBins = unique([ilims(1), i2]);  % Might be distinct or same
    ichunk = [initialBins, (i2+dchunk):dchunk:lastInd];  % Chunk times
    if ichunk(end)~=lastInd
        ichunk(end+1) = lastInd;
    end
    nchunks = length(ichunk)-1;
    % reshape
    ichunk = reshape(repmat(ichunk,2,1),[],1);
    ichunk = reshape(ichunk(2:end-1)',[],nchunks)';
    
    % ==== jchunk: chunks in coarse time ====================================
    jchunk = [floor(ichunk(:,1)*dt)-1, ceil(ichunk(:,2)*dt)+1];
    jchunk(1,1) = max(jchunk(1,1),0);
    jchunk(nchunks,2) = min(OPTprs.slen,jchunk(nchunks,2));

end

% ===== insert params into OPTprs structure =============================
OPTprs.nchunks = nchunks;
OPTprs.ichunk = ichunk;  % fine time bins
OPTprs.jchunk = jchunk;  % coarse time bins

% time interpolation matrix for stimulus
OPTprs.MMintrp = makeInterpMatrix2(max(diff(jchunk')),ndt);

