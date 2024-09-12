function mapReceptiveFields(folder)

%% Parameters
% for correcting baseline drifts of calcium traces at start of experiments
win_decay = 20; % in s, window to test whether baseline is higher than normal
thresh_decay = 1.5; % in std, threshold for drift
win_correct = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
% kernels and RFs simultaneously
lambdas = logspace(-4, 0, 5);
rf_timeLimits = [0.2 0.4];
crossFolds = 10;

%% Loop through all datasets
% Note: for all SC neurons recorded with 2P, recording was in right
% hemisphere and left eye was on video
fprintf('Fit RFs:\n')
subjects = dir(fullfile(folder.data));
subjects = subjects(~startsWith({subjects.name},'.'));

for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        f = fullfile(folder.data, name, date, '001');

        % ignore datasets for which receptive fields were already fit, and
        % datasets for which white noise or circle data do not exist
        if isfile(fullfile(f, "_ss_rf.pValues.npy")) || ...
                ~isfile(fullfile(f, "_ss_sparseNoise.times.npy")) || ...
                ~isfile(fullfile(f, "_ss_circles.intervals.npy"))
            continue
        end

        % load data
        caData = io.getCalciumData(f);
        if isfile(fullfile(f, "_ss_sparseNoise.times.npy"))
            stimData = io.getVisNoiseInfo(f);
        elseif isfile(fullfile(f, "_ss_circles.intervals.npy"))
            stimData = io.getCircleData(f);
        end
        t_stim = stimData.times;
        tBin_stim = median(diff(t_stim));
        t_rf = (floor(rf_timeLimits(1) / tBin_stim) : ...
            ceil(rf_timeLimits(2) / tBin_stim)) .* tBin_stim;

        %% Prepare calcium traces
        % interpolate calcium traces to align all to same time
        t_ind = caData.time > t_stim(1) - 10 & ...
            caData.time < t_stim(end) + 10;
        caTraces = caData.traces(t_ind,:);
        t_ca = caData.time(t_ind);
        [caTraces, t_ca] = traces.alignSampling(caTraces, t_ca, ...
            caData.planes, caData.delays);

        % remove strong baseline decay at start of experiment in cells that
        % show it
        caTraces = traces.removeDecay(caTraces, t_ca, win_decay, ...
            win_correct, thresh_decay);

        % if white noise data available, use it to map RFs
        % otherwise, if circle data available, use that to map RFs
        if isfile(fullfile(f, "_ss_sparseNoise.times.npy"))
            stimFrames = stimData.frames(stimData.stimOrder,:,:);
            
            % map RF
            % rFields: [rows x columns x t_rf x ON/OFF x units]
            [rFields, ev] = ...
                whiteNoise.getReceptiveField(caTraces, t_ca, ...
                stimFrames, t_stim, t_rf./tBin_stim, ...
                lambdas, crossFolds);

            v = squeeze(mean(ev,3)); % [neuron x lamStim]
                                     % average EV across cross-folds
            [maxEV, maxStimLam] = max(v,[],2);
            bestLambdas = lambdas(maxStimLam)';

            % test signficance of each RF (note: resulting ev are not
            % cross-validated, while maxEV are)
            [ev, ev_shift] = ...
                whiteNoise.receptiveFieldShiftTest( ...
                caTraces, t_ca, stimFrames, t_stim, ...
                t_rf./tBin_stim, rFields, bestLambdas, 500);
            pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
            pvals(isnan(ev)) = NaN;

            % save results
            results.maps = rFields;
            results.explVars = maxEV;
            results.lambdas = bestLambdas;
            results.pValues = pvals;
            results.timestamps = t_rf;
            results.edges = stimData.edges;
            io.writeNoiseRFmapResults(results, f);
        
        elseif isfile(fullfile(f, "_ss_circles.intervals.npy"))
            % create stimulus matrix [t_stim x rows x columns x circle_diameter]
            [stimMatrix, x, y, diameters] = circles.getStimMatrix(t_stim, ...
                stimData.xyPos, stimData.diameter, stimData.isWhite);
            % create toeplitz matrix
            [toeplitz, t_toeplitz] = whiteNoise.makeStimToeplitz( ...
                stimMatrix, t_stim, t_rf/tBin_stim);

            % temporal RF for each circle position and size; one RF per eye
            % position bin
            % rFields: [X x Y x diameter x t_rf x ON/OFF x units]
            [rFields, ev] = circles.getReceptiveField(caTraces, t_ca, ...
                toeplitz, t_toeplitz, ...
                [length(x) length(y) length(diameters) length(t_rf)], ...
                lambdas, crossFolds);

            v = squeeze(mean(ev,3)); % [neuron x lamStim]
                                     % average EV across cross-folds
            [maxEV, maxStimLam] = max(v,[],2);
            bestLambdas = lambdas(maxStimLam)';

            % test signficance of each RF (note: resulting ev are not
            % cross-validated, while maxEV are)
            [ev, ev_shift] = ...
                circles.receptiveFieldShiftTest(caTraces, t_ca, ...
                toeplitz, t_toeplitz, rFields, bestLambdas, 500);
            pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
            pvals(isnan(ev)) = NaN;
        end

    end
end