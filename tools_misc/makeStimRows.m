function S= makeStimRows(Stim, n, flag);
%  S = makeStimRows(Stim, n, flag);
%
%  Converts raw stimulus to a matrix where each row is loaded with the full
%  space-time stimulus at a particular moment in time.  The resulting
%  matrix has length equal to the number of time points and width equal to
%  the (number of spatial dimensions) x (kernel size n).
%
%  Inputs: 
%   Stim = stimulus, first dimension is time, other dimensions are spatial
%          dimensions
%   n = size of temporal kernel; number of time samples to include in each
%       row of stimulus matrix.
%   flag (optional)
%        = 0, default behavior: padded w/ zeros at beginning so length of
%          output matrix matches size of Stim
%        = 1, no padding with zeros: length of S is length(Stim)-n+1.
%        = vector of indices, (e.g. indices of spikes times).  Return
%        a matrix containting only the spiking stimuli
%
%  Output: S = matrix where each row is the size of the linear kernel
%
%  Last updated:  6/30/2005, J. Pillow

% parse inputs
if nargin < 3
    flag = 0;
    if nargin < 2
        global n
        if isempty(n)
            error('ERROR -- makeStimRows:  n is undefined');
        end
    end
end

sz = size(Stim);
n2 = prod(sz(2:end));  % total dimensionality in spatial dimensions

% If necessary, convert Stim to a 2D matrix
if (n2 > sz(2));       % reshape to matrix if necessary
    Stim = reshape(Stim, sz(1), n2);
end

if flag == 0  % Compute with zero-padding. ----------------------------------
    S = zeros(sz(1), n2*n);
    for j=1:n2
        S(:,(n*(j-1)+1):(n*j)) = ...
            fliplr(toeplitz(Stim(:, j), [Stim(1,j) zeros(1,n-1)]));
    end
    
    
elseif (length(flag) == 1)  % compute only for stimuli at times >= n ----------
    S = zeros(sz(1)-n+1,n2*n);
    for j=1:n2
        S(:,(n*(j-1)+1):(n*j)) = ...
            fliplr(toeplitz(Stim(n:sz(1), j), Stim(n:-1:1,j)));
    end
    
    
else % compute for spiking stimuli --------------------------------------------
    if (min(flag) < 1) | (max(flag) > sz(1))
        error('makeStimRows:  3rd arg should be spike indices (vals are too high or too low): %d %d', ...
            min(flag), max(flag));
    end
    S = zeros(length(flag), n2*n);
    % Do for spinds < n
    nsp1 = length(flag(flag<n));
    for j = 1:nsp1
        S(j,:) = reshape([zeros(n-flag(j), n2); Stim(1:flag(j),:)], 1, n2*n);
    end
    for j = nsp1 +1:length(flag)
        S(j,:) = reshape(Stim(flag(j)-n+1:flag(j),:), 1, n2*n);
    end
end
