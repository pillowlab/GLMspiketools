function Y = spikefilt(sps, H)
% Y = spikefilt(sps, H, twin)
%
% Causally filters spike train with spike indices spikeInds with filter matrix H.
%
% Inputs:
%   sps [tlen x 1] - spike train vector
%       H  [n x m] - spike-history matrix (each column is a filter)
%
% Output: 
%    Y [tlen x m] - each column is the spike train filtered with the
%                   corresponding column of H
%
% Note: there is a faster 'mex' version of this function called spikefilt_mex,
% available in version 'v1'.  It was removed in Jan 2016 to allow for
% mex-free code release.


[hlen,hwid] = size(H);

% Do convolution and remove extra bins
Y = conv2(sps,H,'full');
Y = [zeros(1,hwid); Y(1:end-hlen,:)];
