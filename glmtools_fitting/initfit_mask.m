function [iiSpk, iiLi] = initfit_mask(mask,dt,rlen)
% [iiSpk, iiLi] = initfit_mask(mask,rlen)
%
% Extracts spike times and bin indices for computing likelihood, given
% intervals contained in mask (an n x 2 matrix)

global SPNDS OPTprs

mask = round(mask./dt);

% Check that mask doesn't extend beyond stimulus window
if max(mask(:)) > rlen
    warning('GLM:mask', 'Mask is too long for data segment: truncating...');
    rowmax = max(mask,[],2);
    maskkeep = (rowmax <= rlen);
    mask = mask(maskkeep,:);
end

% Generate mask
if isempty(mask)  % No mask
    iiLi = [1:rlen];
    iiSpk = SPNDS;
else  % Compute mask for time bins and spikes to use
    masklens = diff(mask,1,2);
    iiLi = zeros(sum(masklens),1);
    icum = 0;
    for j = 1:size(mask,1)
        iiLi(icum+1:icum+masklens(j)) = (mask(j,1)+1:mask(j,2))';
        icum = icum+masklens(j);
    end
    iiSpk = intersect(SPNDS,iiLi);
end

% Put params into global optimization structure OPTprs
OPTprs.mask = mask;
OPTprs.iiSpk = iiSpk;
OPTprs.iiLi = iiLi;
