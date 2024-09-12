function fitRFFromCircles(folder)

%% Parameters
% for binning of horizontal eye position
numEyeBins = 10;

% for correcting baseline drifts of calcium traces at start of experiments
win_decay = 20; % in s, window to test whether baseline is higher than normal
thresh_decay = 1.5; % in std, threshold for drift
win_correct = 150; % in s, window to fit exponential

% to fit RFs from circle protocol
rf_timeLimits = [0 0.5];

%% Loop through all datasets
% Note: for all SC neurons recorded with 2P, recording was in right
% hemisphere and left eye was on video
fprintf('Fit RFs across pupil positions:\nDatasets:\n')

subjects = dir(folder.data);
subjects = subjects(~startsWith({subjects.name},'.'));
for subj = 1 %1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        f = fullfile(folder.data, name, date, '001');
        % only consider datasets where circles were presented
        if ~isfile(fullfile(f, '_ss_circles.intervals.npy'))
            continue
        end
        
        % load data
        caData = io.getCalciumData(f);
        pupilData = io.getPupilData(f);
        stimData = io.getCircleData(f);

        % interpolate pupil data to match stimulus times
        t_stim = stimData.times;
        t_eye = pupilData.time;
        ind = t_eye > t_stim(1)-5 & ...
            t_eye < t_stim(end)+5 & [1;diff(t_eye)] > 0;
        t_eye = t_eye(ind);
        xFaster = ceil(median(diff(t_stim), "omitnan") / ...
            median(diff(t_eye), "omitnan"));
        eyePosX = smooth(pupilData.pos(ind,1), xFaster);
        eyePosX = interp1(t_eye, eyePosX, t_stim);

        % get median eye position for each bin
        [~,~,bin] = histcounts(eyePosX, ...
            prctile(eyePosX, linspace(0,100,numEyeBins+1)));
        eyePosPerBin = NaN(1, numEyeBins);
        for b = 1:numEyeBins
            eyePosPerBin(b) = median(eyePosX(bin == b), "omitnan");
        end

        %% Prepare calcium traces
        % interpolate calcium traces to align all to same time
        t_ind = caData.time > t_stim(1) - 10 & caData.time < t_stim(end) + 10;
        caTraces = caData.traces(t_ind,:);
        t_ca = caData.time(t_ind);
        [caTraces, t_ca] = traces.alignSampling(caTraces, t_ca, ...
            caData.planes, caData.delays);

        % remove strong baseline decay at start of experiment in cells that
        % show it
        caTraces = traces.removeDecay(caTraces, t_ca, win_decay, ...
            win_correct, thresh_decay);

        %% Map shifted RFs
        tBin_stim = median(diff(t_stim));
        t_rf = (floor(rf_timeLimits(1)/tBin_stim) : ...
            ceil(rf_timeLimits(2)/tBin_stim)) .* tBin_stim;

        % create stimulus matrix [t_stim x rows x columns x circle_diameter]
        [stimMatrix, x, y, diameters] = circles.getStimMatrix(t_stim, ...
            stimData.xyPos, stimData.diameter, stimData.isWhite);

        % create toeplitz matrix
        [toeplitz, t_toeplitz] = whiteNoise.makeStimToeplitz(stimMatrix, t_stim, ...
            t_rf/tBin_stim);

        % temporal RF for each circle position and size; one RF per eye
        % position bin
        % [X x Y x diameter x t_rf x ON/OFF x units]
        rfs = circles.getReceptiveField(caTraces, t_ca, toeplitz, ...
            t_toeplitz, ...
            [length(x) length(y) length(diameters) length(t_rf)], ...
            lambda, ignoreTimes);
        % positive pixel in OFF subfield means: driven by black
        rfs(:,:,:,:,2,:) = -rfs(:,:,:,:,2,:);

        for iUnit = 1:size(caTraces,2)
            % find best time (for white and black separately)
            max(rfs,[],[1 2 3])

            % find best diameter (for white and black separately)

            % find center of mean of ON- and OFF-subfields
        end

        %% Linear regression
        % linear regression relating pupil position in pixels to RF centre
        % in visual degrees &
        % determine median horizontal RF positions for newly fitted RFs,
        % only use units where all of the shifted x-positions are within
        % stimulus

        %% Save results

        %% Make plots
    end
end