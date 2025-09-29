function [zTraces, stim, validTimes, validUnits, validStim, ...
    lamStim, lamMatrix_stim] = ...
    prepareDataForRFFit(caTraces, t_ca, toeplitz, t_toeplitz, ...
    stimSize, rfBins, lambdas, ignoreStimTimes)

% get neural response
tBin_ca = median(diff(t_ca));
tBin_stim = median(diff(t_toeplitz));
numBins = round(tBin_stim / tBin_ca);
caTraces = smoothdata(caTraces, 1, 'movmean', numBins, 'omitnan');
% resample neural response at stimulus times
caTraces = interp1(t_ca, caTraces, t_toeplitz);
% z-score neural response
zTraces = (caTraces - mean(caTraces,1,'omitnan')) ./ ...
    std(caTraces,0,1,'omitnan');

% clean up neural traces (delete times where all traces are NaN; if NaN 
% values < 10% in a neuron, exchange NaNs for 0; skip neurons that have only 
% NaN values)
validTimes = ~all(isnan(zTraces),2);
toeplitz(~validTimes,:) = [];
zTraces(~validTimes,:) = [];
ignoreStimTimes(~validTimes) = [];
% if NaN values < 10% in a neuron, exchange NaNs for 0
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.1;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
% skip neurons that have only NaN values
validUnits = ~any(isnan(zTraces),1)';

if any(ignoreStimTimes)
    % set stimuli that occured at times ignoreStimTimes to zero; Note: these
    % are not the times of the neural response but the times of the stimulus
    % occurrence; this is used, for example, to only consider stimuli that were
    % presented during specific pupil positions
    stimPars = prod(stimSize);
    for b = 1:length(rfBins)
        ignore = [true(rfBins(b),1); ignoreStimTimes(1:end-rfBins(b))];
        toeplitz(ignore, (b-1)*stimPars + (1:stimPars)) = 0;
    end
end
% delete stimulus frames that will be ignored
indTime = all(toeplitz==0,2);
indVal = find(validTimes);
validTimes(indVal(indTime)) = false;
zTraces(indTime,:) = [];
toeplitz(indTime,:) = [];

% duplicate stimulus matrix to predict ON part (1st half) and OFF
% part (2nd half)
s = toeplitz;
s(toeplitz < 0) = 0;
stim = s;
s = toeplitz;
s(toeplitz > 0) = 0;
stim = [stim, s];
validStim = ~all(stim==0,1);
% normalise each column of stimulus matrix
stim = (stim - mean(stim(:),'omitnan')) ./ std(stim(:),'omitnan');

% scale lambdas according to number of samples and number of predictors
lamStim = sqrt(lambdas .* sum(validTimes) .* size(stim,2));
% construct spatiotemporal smoothing lambda matrix
lamMatrix_stim = krnl.makeLambdaMatrix([stimSize, length(rfBins)], ...
    ones(1,length(stimSize)+1));
lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);