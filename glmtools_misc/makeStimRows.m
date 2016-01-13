function S = makeStimRows(Stim, nkt, flag)
% S = makeStimRows(Stim, nkt, flag);
%
% Converts spatio-temporal stimulus to a design matrix, where each row of
% the design matrix contains all space-time stimulus elements within a time
% window of "nkt" time bins. 
%
% INPUT: 
%  Stim [N x M] - stimulus, first dimension is time, others are space
%   nkt [1 x 1] - number of time bins to include in temporal kernel
%  flag (optional)
%         'same' - padded w/ zeros to make S the same size (default)
%         'valid' - no padding with zeros: length(S) = length(Stim)-nkt+1.
%          vector of indices - returns only stimuli corresponding to
%                                these time bins 
%
% OUTPUT: 
%  S [N x M*nkt] or  [Nsp x M*nkt ] - design matrix 
%   
%  $Id$

% parse inputs
if nargin < 3
    flag = 'same';
end

[n1,n2] = size(Stim);

if ischar(flag) % --- Form for all stimuli (with zero-padding of first bins) -----

    if strcmp(flag,'same')
        S = zeros(n1, n2*nkt);
        Stim = [zeros(nkt-1,n2); Stim];
        n1 = n1+nkt-1;
    elseif strcmp(flag,'valid')
        S = zeros(n1-nkt+1,n2*nkt);
    else
        error('Unknown flag');
    end
    
    % Form design matrix
    for j=1:n2
        % Loop over columns of original stimulus
        S(:,(nkt*(j-1)+1):(nkt*j)) = ...
            fliplr(toeplitz(Stim(nkt:n1, j), Stim(nkt:-1:1,j)));
    end
    
    
else % ----  Form for spiking stimuli only ----------------------------------------
    inds = flag; % indices to keep
    if (min(inds) < 1) || (max(inds) > n1)
        error('Spike indices outside allowed range');
    end
    S = zeros(length(inds), n2*nkt);

    % Do for spinds < nkt
    nsp1 = length(inds(inds<nkt));
    for j = 1:nsp1
        S(j,:) = reshape([zeros(nkt-inds(j), n2); Stim(1:inds(j),:)], 1, n2*nkt);
    end
    
    % Do for remaining indices
    for j = nsp1 +1:length(inds)
        S(j,:) = reshape(Stim(inds(j)-nkt+1:inds(j),:), 1, n2*nkt);
    end
end
