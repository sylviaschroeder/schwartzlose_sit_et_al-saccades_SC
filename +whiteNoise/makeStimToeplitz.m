function [stim, time, stimFrames, stimBin] = ...
    makeStimToeplitz(stimFrames, stimTimes, RFtimesInFrames)

%MAKESTIMTOEPLITZ Generates toeplitz matrix for whitenoise stimuli.
%   [stim, time, stimFrames, stimBin] = ...
%    MAKESTIMTOEPLITZ(stimFrames, stimTimes, RFtimesInFrames) generates the
%    toeplitz matrix for the given stimulus frames.
%
%   stim                [time x pixels]; toeplitz matrix
%   time                [time x 1]; time of stimulus frames (no gaps)
%   stimFrames          [time x rows x cols]; stimulus frames, gaps filled
%                       with zeros
%   stimBin             double; duration of each stimulus frame
%
%   stimFrames          [time x rows x cols]; noise stimulus
%   stimTimes           [time x 1]; times of stimulus frames
%   RFtimesInFrames     [1 x RFframes]; frames of receptive field relative
%                       to stimulus frames

% generate toeplitz matrix for stimuli: [time x pixels]
% each row holds all pixels at current and previous time points:
% [[all pixels at t=0], [all pixels at t=-1], ...]
% each column is time series of that particular pixel

% find time gaps in stimulus presentation (usually when same visual noise
% stimulus was repeated several times)
stimBin = median(diff(stimTimes));
indGap = find(diff(stimTimes) > 2 * stimBin);
time = stimTimes;
% fill gaps with zeros in stimulus matrix
for g = 1:length(indGap)
    add = round(diff(stimTimes(indGap(g) + [0 1])) ./ stimBin);
    stimFrames = [stimFrames(1:indGap(g),:,:); ...
        zeros(add, size(stimFrames,2), size(stimFrames,3)); ...
        stimFrames(indGap(g)+1:end,:,:)];
    time = [time(1:indGap(g)); ...
        time(indGap(g)) + (1:add)' .* stimBin; ...
        time(indGap(g)+1:end)];
end
% reshape stimulus frames to [time x px]; this represents a single
% "stimulus block", i.e. the pixels to estimate a single time point of the
% receptive field
stim = reshape(stimFrames, size(stimFrames,1), []);
% now concatinate time shifted stimulus blocks; for each time point there
% is a stimulus block for lag=0, another for lag=-1, another for lag=-2,...
st = [];
for t = 1:length(RFtimesInFrames)
    st = [st, ...
        [zeros(max(0,RFtimesInFrames(1)-1+t), size(stim,2)); ...
        stim(max(1,2-RFtimesInFrames(1)-t) : end-RFtimesInFrames(1)-t+1, :)]];
end
stim = st;