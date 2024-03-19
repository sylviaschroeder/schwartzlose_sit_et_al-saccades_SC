%% Parameters
% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;
% for binning of horizontal eye position
numEyeBins = 10;

%% Fit RF to each eye position for cells with good overall RFs
% Note: for all SC neurons recorded with 2P, recording was in right
% hemisphere and left eye was on video
subjects = dir(fullfile(folderData, 'SS*'));
for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folderData, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        folder = fullfile(folderData, name, date, '001');

        if ~isfile(fullfile(folder, '_ss_sparseNoise.times.npy'))
            continue
        end
        
        % load data
        caData = io.getCalciumData(folder);
        pupilData = io.getPupilData(folder);
        noiseData = io.getVisNoiseInfo(folder);
        rfData = io.getRFData(folder);

        % interpolate pupil data to match stimulus times

        % smooth pupil data! (unless already done)

        t_stim = noiseData.times;
        eyePos = interp1(pupilData.time, pupilData.pos(:,1), t_stim);

        % percentiles of horizontal eye positions
        % Note: psition from nasal (smaller) to temporal (larger)
        eyePosQuart = prctile(pupilData.pos(:,1), [25 50 75]);
        eyePosNorm = (pupilData.pos(:,1) - eyePosQuart(2)) ./ ...
            diff(eyePosQuart([1 3]));
        eyeStart = find(pupilData.time >= noiseData.interval(1), 1);
        eyeEnd = find(pupilData.time > noiseData.interval(2), 1) - 1;
        [~,eyeBinEdges,bin] = histcounts(eyePosNorm(eyeStart:eyeEnd), ...
            prctile(eyePosNorm(eyeStart:eyeEnd), 0:numEyeBins:100));

        validRF = find(rfData.pValues < minPVal & ...
            rfData.explVars > minExplainedVarianceStim & ...
            rfData.lambdas < maxLambda);
        for k = 1:length(validRF)
            iCell = validRF(k);

        end
    end
end