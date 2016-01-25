function Y = spikefilt(spikeInds, H, twin)
% Y = spikeilft(spikeInds, H, twin)
%
% Caussally filters spike train with spike indices spikeInds with filter matrix H.
%
% This is a slower m-file version of the mex file spikefilt_mex.
% Users with difficulty compiling mex can use this function instead.

[hlen,hwid] = size(H);

% Create spike train vector
spnds = in(spikeInds,twin)-twin(1)+1;
sps = zeros(diff(twin)+1,1);
sps(spnds) = 1;

% Do convolution and remove extra bins
Y = conv2(sps,H,'full');
Y = [zeros(1,hwid); Y(1:end-hlen,:)];
