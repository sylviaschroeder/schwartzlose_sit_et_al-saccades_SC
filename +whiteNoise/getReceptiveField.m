function [receptiveFields, explainedVariance, predictions] = ...
    getReceptiveField(caTraces, t_ca, ...
    toeplitz, t_toeplitz, stimSize, rfBins, ...
    lambdas, crossFolds, ignoreStimTimes)

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

if nargin < 9
    ignoreStimTimes = false(size(t_toeplitz));
end

[zTraces, stim, validTimes, validUnits, validStim, ...
    lamStim, lamMatrix_stim] = ...
    whiteNoise.prepareDataForRFFit(caTraces, t_ca, toeplitz, t_toeplitz, ...
    stimSize, rfBins, lambdas, ignoreStimTimes);

% for cross-validation, choose every n-th sample for one "fold" (where n is
% number of folds); don't choose chunks of successive samples
nPerFold = ceil(size(stim,1) / crossFolds);
indPerFold = reshape(1:(crossFolds*nPerFold), crossFolds, [])';
indPerFold(indPerFold > size(stim,1)) = NaN;

explainedVariance = NaN(size(caTraces,2), length(lamStim), crossFolds);
preds = NaN(nPerFold, crossFolds, size(caTraces,2), length(lamStim));

% get variances explained
for fold = 1:crossFolds
    ind = indPerFold(:,fold);
    ind(isnan(ind)) = [];
    j = true(size(zTraces,1),1);
    j(ind) = false;
    
    if crossFolds > 1
        y_train = gpuArray(padarray(zTraces(j,validUnits), ...
            size(lamMatrix_stim,1), 'post'));
        y_mean = mean(zTraces(j, validUnits),1);
        x_train = stim(j,:);
    else
        y_train = gpuArray(padarray(zTraces(~j,validUnits), ...
            size(lamMatrix_stim,1), 'post'));
        y_mean = mean(zTraces(~j, validUnits),1);
        x_train = stim(~j,:);
    end
    y_test = zTraces(~j,validUnits);
    x_test = stim(~j,:);

    for lamS = 1:length(lamStim)
        lms = lamMatrix_stim .* lamStim(lamS);
        
        A = gpuArray([x_train; lms]);

        B = gather(A \ y_train);
        pred = x_test * B; % get prediction
        preds(1:size(pred,1), fold, validUnits, lamS) = pred;
        explainedVariance(validUnits, lamS, fold) = 1 - ...
            sum((y_test - pred) .^ 2,1) ./ ...
            sum((y_test - y_mean) .^2, 1);
    end
end

receptiveFields = NaN(size(stim,2), size(caTraces,2));
if length(lamStim) > 1 || crossFolds > 1
    % determine RFs using all data and optimal lambdas
    [~, bestStimLams] = max(mean(explainedVariance, 3), [], 2);
    for lamS = 1:length(lamStim)
        ind = bestStimLams == lamS & validUnits;
        if sum(ind) == 0
            continue
        end
        A = [stim; lamMatrix_stim .* lamStim(lamS)];
        tr = padarray(zTraces(:,ind), size(lamMatrix_stim,1), 'post');

        B = gather(gpuArray(A) \ gpuArray(tr));
        receptiveFields(:,ind) = B; % get RF kernel
        preds(:,:,ind,1) = preds(:,:,ind,lamS);
    end
else
    receptiveFields(:,validUnits) = B;
end

receptiveFields(~validStim,:) = NaN;
receptiveFields = reshape(receptiveFields, ...
    [stimSize, length(rfBins), 2, size(caTraces,2)]);

preds = reshape(preds(:,:,:,1), [], size(caTraces,2));
preds = preds(1:size(stim,1),:);
predictions = NaN(length(t_toeplitz), size(caTraces,2));
predictions(validTimes,:) = preds;