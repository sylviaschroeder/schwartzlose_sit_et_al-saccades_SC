function [rfs, explainedVariance, predictions] = getReceptiveField( ...
    caTraces, t_ca, toeplitz, t_toeplitz, ...
    stimSize, lambdas, crossFolds, ignoreTimes)

if nargin < 8
    ignoreTimes = false(size(t_toeplitz));
end

% resample neural responses at stimulus times
tBin_ca = median(diff(t_ca));
tBin_stim = median(diff(t_toeplitz));
xTimes = ceil(tBin_stim / tBin_ca);
caTraces = smoothdata(caTraces, 1, 'movmean', xTimes, 'omitnan');
caTraces = interp1(t_ca, caTraces, t_toeplitz);

% z-score neural responses
zTraces = (caTraces - mean(caTraces,1,"omitnan")) ./ ...
    std(caTraces,0,1,"omitnan");

% clean up neural traces (delete times where all traces are NaN; if NaN 
% values < 5% in a neuron, exchange NaNs for 0; skip neurons that have only 
% NaN values)
validTimes = ~all(isnan(zTraces),2);
toeplitz(~validTimes,:) = [];
zTraces(~validTimes,:) = [];
ignoreTimes(~validTimes) = [];
ind = any(isnan(zTraces),1) & sum(isnan(zTraces),1)/size(zTraces,1) <= 0.1;
if sum(ind) > 0
    zTraces(:,ind) = fillmissing(zTraces(:,ind),'constant',0);
end
validUnits = ~any(isnan(zTraces),1)';

% duplicate stimulus matrix to predict ON part (1st half) and OFF
% part (2nd half)
s = toeplitz;
s(toeplitz < 0) = 0;
stim = s;
s = toeplitz;
s(toeplitz > 0) = 0;
stim = [stim, s];
% normalise each column of stimulus matrix
stim = (stim - mean(stim(:),'omitnan')) ./ std(stim(:),'omitnan');

% for stimulus frames that will be ignored, set neural traces and stimulus
% to zero
zTraces(ignoreTimes,:) = 0;
stim(ignoreTimes,:) = 0;

% scale lambdas and construct smoothing matrix
lamStim = sqrt(lambdas .* size(stim,1) .* size(stim,2));
lamMatrix_stim = krnl.makeLambdaMatrix(stimSize, [1 1 1 1]);
% duplicate to cover white and black circles
lamMatrix_stim = blkdiag(lamMatrix_stim, lamMatrix_stim);

% for cross-validation, choose every n-th sample for one "fold" (where n is
% number of folds); don't choose chunks of successive samples
nPerFold = ceil(size(stim,1) / crossFolds);
indPerFold = reshape(1:(crossFolds*nPerFold), crossFolds, [])';
indPerFold(indPerFold > size(stim,1)) = NaN;

explainedVariance = NaN(size(caTraces,2), length(lamStim), crossFolds);
preds = NaN(nPerFold, crossFolds, size(caTraces,2), length(lamStim));

% get variances explained
for fold = 1:crossFolds
    % ind = (1:nPerFold) + (fold-1)*nPerFold;
    % ind(ind > size(zTraces,1)) = [];
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
        preds(1:length(pred), fold, validUnits, lamS) = pred;
        explainedVariance(validUnits, lamS, fold) = 1 - ...
            sum((y_test - pred) .^ 2,1) ./ ...
            sum((y_test - y_mean) .^2, 1);
    end
end

if length(lamStim) > 1 || crossFolds > 1
    % determine RFs using all data and optimal lambdas
    rfs = NaN(size(stim,2), size(caTraces,2));

    [~, bestStimLams] = max(mean(explainedVariance, 3), [], 2);
    for lamS = 1:length(lamStim)
        ind = bestStimLams == lamS & validUnits;
        if sum(ind) == 0
            continue
        end
        A = [stim; lamMatrix_stim .* lamStim(lamS)];
        tr = padarray(zTraces(:,ind), size(lamMatrix_stim,1), 'post');

        B = gather(gpuArray(A) \ gpuArray(tr));
        rfs(:,ind) = B; % get RF kernel
        preds(:,:,ind,1) = preds(:,:,ind,lamS);
    end
else
    rfs = B;
end

rfs = reshape(rfs, [stimSize 2 size(caTraces,2)]);

preds = reshape(preds(:,:,:,1), [], size(caTraces,2));
preds = preds(1:size(stim,1),:);
predictions = NaN(length(t_toeplitz), size(caTraces,2));
predictions(validTimes,:) = preds;