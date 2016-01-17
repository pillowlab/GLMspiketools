function bmask = initfit_mask(mask,dtSp,rlen)
% bmask = initfit_mask(mask,dtSp,rlen)
%
% Compute binary mask for computing log-likelihood
%
% Input:
% ------
%   mask [n x 2]: list of intervals to use for computing log-likelihood
%                   (ignore time bins outside these intervals)
%   dtSp [1 x 1]: bin size for spike trains and conditional intensity
%   rlen [1 x 1]: length of total conditional intensity vector
%
% Output: 
% -------  
%  bmask [rlen x 1]: binary vector where (1 = use for LL), (0 = ignore).

if isempty(mask)  % No mask
    bmask = true(rlen,1); % keep all indices

else 
    
    % convert mask to discrete time bins
    mask = round(mask/dtSp);

    % Check that mask doesn't extend beyond stimulus window
    if max(mask(:)) > rlen
        warning('GLM:mask', 'Mask is too long for data segment: truncating...');
        mask(:,2) = max(mask(:,2),rlen); % set maximum mask value to stim length
        mask(diff(mask')'<0,:) = []; % remove rows with negative intervals
    end
    
    % Generate mask
    bmask = false(rlen,1);
    for j = 1:size(mask,1)
        bmask(mask(j,1):mask(j,2)) = true;
    end
end
