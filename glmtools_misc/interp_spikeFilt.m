function ih1 = interp_spikeFilt(ih0,iht,dt)
% ih1 = interp_spikeFilt(ih0,iht,dt)
%
% If necessary, interpolate post-spike filter to correct time binning and
% cut off any bins <= 0 (i.e., the spike time itself).
%
% Inputs:  ih0 = post-spike waveform (or matrix of post-spike waveforms)
%          iht = time binning of ih0
%          dt = desired time binning
%
% Outputs ih1 = interpolated, possibly-truncated version of ih0

if ~isempty(ih0)
    htest = round(10*diff(iht)./dt);  % Check if time bins for h same as dt 
    if all(htest == 10)  
        ih1 = ih0(iht>0,:); % cut off times <= 0
    else % Interpololate for binning at dt
        iht1 = [dt:dt:max(iht)]';  % time points for sampling
        ih1 = interp1(iht, ih0, repmat(iht1,1,size(ih0,2)), 'linear', 0);
    end

else % ih0 is empty!
    ih1 = zeros(1,size(ih0,2));
    if ~isempty(iht)
        fprintf('Warning: no post-spike kernel supplied but iht nonempty!!\n');
    end
end
