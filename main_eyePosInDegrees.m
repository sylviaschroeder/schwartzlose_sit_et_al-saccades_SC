%% Parameters
% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;

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

        % percentiles of horizontal eye positions
        % Note: psition from nasal (smaller) to temporal (larger)
        eyePosPerc = prctile(pupilData.pos(:,1), 5:10:95) - ...
            median(pupilData.pos(:,1), 'omitnan');

        validRF = find(rfData.pValues < minPVal & ...
            rfData.explVars > minExplainedVarianceStim & ...
            rfData.lambdas < maxLambda);
        for k = 1:length(validRF)
            iCell = validRF(k);

        end
    end
end