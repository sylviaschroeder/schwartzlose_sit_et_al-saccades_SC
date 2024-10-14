function [receptiveFields, explainedVariance, predictions, time] = ...
    getReceptiveField(traces, traceTimes, ...
    stimFrames, stimTimes, RFtimesInFrames, ...
    lambdas, crossFolds, ignoreFrames)

%GETRECEPTIVEFIELD Returns spatiotemporal receptive field.
%   [receptiveFields, explainedVariance, predictions, time] = ...
%    GETRECEPTIVEFIELD(traces, traceTimes, ...
%    stimFrames, stimTimes, RFtimesInFrames, ...
%    lambdas, crossFolds) calculates the linear RF of the neuron.
%
%   receptiveFields     [rows x cols x RFframes x RFtype x neuron]
%                       containing linear regression solution for x in Ax=B
%                       where A is stimulus [rows x cols x time] and B is
%                       calcium response, for each neuron and stimulus 
%                       model; ridge regression is performed on all data
%                       using the optimal lambda value found with
%                       cross-validation
%   explainedVariance   [neuron x lambdaStim x crossFold], each entry:
%                       explained variance for fitted RF for
%                       each neuron, lambda, and cross val. fold
%   predictions         [t x neuron], each column contains
%                       prediction based on RF for test 
%                       responses of specific neuron (using optimal
%                       lambda)
%   time                [t x 1]; time points for predictions
%
%   traces              [trTime x neuron]; calcium traces of neurons
%   traceTimes          [trTime x 1]; sample times of calcium traces
%   stimFrames          [time x rows x cols]; noise stimulus
%   stimTimes           [time x 1]; times of stimulus frames
%   RFtimesInFrames     [1 x RFframes]; frames of receptive field relative
%                       to stimulus frames
%   lambdas             [1 x lambda]; values of lambda
%   crossFolds          ind; number of cross val. folds

% generate toeplitz matrix for stimulus
[stim, time, stimBin] = ...
    whiteNoise.makeStimToeplitz(stimFrames, stimTimes, RFtimesInFrames);

% get neural response
traceBin = median(diff(traceTimes));
numBins = round(stimBin / traceBin);
traces = smoothdata(traces, 1, 'movmean', numBins, 'omitnan');
% resample neural response at stimulus times
traces = interp1(traceTimes, traces, time);
% z-score neural response
zTraces = (traces - mean(traces,1,'omitnan')) ./ std(traces,0,1,'omitnan');

% delete stim frames for which all neurons have NaN
ind = all(isnan(zTraces),2);
stim(ind,:) = [];
zTraces(ind,:) = [];
time(ind) = [];
% if NaN values < 5% in a neuron, exchange NaNs for 0
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.1;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
% skip neurons that have only NaN values
valid = ~any(isnan(zTraces),1)';

% duplicate stimulus matrix to predict ON part (1st half) and OFF
% part (2nd half)
s = stim;
s(stim < 0) = 0;
stim2 = s;
s = stim;
s(stim > 0) = 0;
stim2 = [stim2, s];
% normalise each column of stimulus matrix
stim2 = (stim2 - mean(stim2(:),'omitnan')) ./ std(stim2(:),'omitnan');

% for stimulus frames that will be ignored, set neural traces and stimulus
% to zero
if nargin > 7
    % find time points to be ignored
    ind = ismember(time, stimTimes(ignoreFrames));
    % set traces and stimulus to zero
    zTraces(ind,:) = 0;
    stim2(ind,:) = 0;
end

if isempty(lambdas)
    lamStim = 0;
    lamMatrix_stim = [];
else
    % scale lambdas according to number of samples and number of predictors
    lamStim = sqrt(lambdas .* size(stim,1) .* size(stim,2));

    % construct spatial smoothing lambda matrix
    lamMatrix_stim = krnl.makeLambdaMatrix([size(stimFrames,2), size(stimFrames,3), ...
        length(RFtimesInFrames)], [1 1 1]);
    lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);
end

% for cross-validation, choose every n-th sample for one "fold" (where n is
% number of folds); don't choose chunks of successive samples
nPerFold = ceil(size(stim,1) / crossFolds);
indPerFold = reshape(1:(crossFolds*nPerFold), crossFolds, [])';
indPerFold(indPerFold > size(stim,1)) = NaN;

explainedVariance = NaN(size(traces,2), length(lamStim), crossFolds);
predictions = NaN(nPerFold, crossFolds, size(traces,2), length(lamStim));

% get variances explained
for fold = 1:crossFolds
    % ind = (1:nPerFold) + (fold-1)*nPerFold;
    % ind(ind > size(zTraces,1)) = [];
    ind = indPerFold(:,fold);
    ind(isnan(ind)) = [];
    j = true(size(zTraces,1),1);
    j(ind) = false;
    
    if crossFolds > 1
        y_train = gpuArray(padarray(zTraces(j,valid), size(lamMatrix_stim,1), 'post'));
        y_mean = mean(zTraces(j, valid),1);
        x_train = stim2(j,:);
    else
        y_train = gpuArray(padarray(zTraces(~j,valid), size(lamMatrix_stim,1), 'post'));
        y_mean = mean(zTraces(~j, valid),1);
        x_train = stim2(~j,:);
    end
    y_test = zTraces(~j,valid);
    x_test = stim2(~j,:);

    for lamS = 1:length(lamStim)
        lms = lamMatrix_stim .* lamStim(lamS);
        
        A = gpuArray([x_train; lms]);

        B = gather(A \ y_train);
        pred = x_test * B; % get prediction
        predictions(1:size(pred,1), fold, valid, lamS) = pred;
        explainedVariance(valid, lamS, fold) = 1 - ...
            sum((y_test - pred) .^ 2,1) ./ ...
            sum((y_test - y_mean) .^2, 1);
    end
end

receptiveFields = NaN(size(stim2,2), size(traces,2));
if length(lamStim) > 1 || crossFolds > 1
    % determine RFs using all data and optimal lambdas
    [~, bestStimLams] = max(mean(explainedVariance, 3), [], 2);
    for lamS = 1:length(lamStim)
        ind = bestStimLams == lamS & valid;
        if sum(ind) == 0
            continue
        end
        A = [stim2; lamMatrix_stim .* lamStim(lamS)];
        tr = padarray(zTraces(:,ind), size(lamMatrix_stim,1), 'post');

        B = gather(gpuArray(A) \ gpuArray(tr));
        receptiveFields(:,ind) = B; % get RF kernel
        predictions(:,:,ind,1) = predictions(:,:,ind,lamS);
    end
else
    receptiveFields(:,valid) = B;
end

receptiveFields = reshape(receptiveFields, size(stimFrames,2), ...
    size(stimFrames,3), length(RFtimesInFrames), 2, size(traces,2));

predictions = reshape(predictions(:,:,:,1), [], size(traces,2));
predictions = predictions(1:length(time),:);