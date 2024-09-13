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

% for evaluation of receptive fields (significance/goodness)
minExplainedVariance = 0.01;
maxPVal = 0.05;
% for plotting RFs
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
titles = {'ON field','OFF field'};

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

        %% Map RFs
        % if white noise data available, use it to map RFs
        % otherwise, if circle data available, use that to map RFs
        if isfile(fullfile(f, "_ss_sparseNoise.times.npy"))
            % load data
            caData = io.getCalciumData(f);
            stimData = io.getVisNoiseInfo(f);
            t_stim = stimData.times;
            tBin_stim = median(diff(t_stim));
            t_rf = (floor(rf_timeLimits(1) / tBin_stim) : ...
                ceil(rf_timeLimits(2) / tBin_stim)) .* tBin_stim;

            % prepare calcium traces
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
            results.maps = permute(rFields, [5 1 2 3 4]);
            results.explVars = maxEV;
            results.lambdas = bestLambdas;
            results.pValues = pvals;
            results.timestamps = t_rf;
            results.edges = stimData.edges;
            io.writeNoiseRFResults(results, f);
        end
        if isfile(fullfile(f, "_ss_circles.intervals.npy"))
            % load data
            caData = io.getCalciumData(f);
            stimData = io.getCircleData(f);
            t_stim = stimData.times;
            tBin_stim = median(diff(t_stim));
            t_rf = (floor(rf_timeLimits(1) / tBin_stim) : ...
                ceil(rf_timeLimits(2) / tBin_stim)) .* tBin_stim;

            % prepare calcium traces
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

            % save results
            results.maps = permute(rFields, [6 1 2 3 4 5]);
            results.explVars = maxEV;
            results.lambdas = bestLambdas;
            results.pValues = pvals;
            results.timestamps = t_rf;
            results.x = x;
            results.y = y;
            results.diameters = diameters;
            io.writeCircleRFResults(results, f);
        end

        %% Plot RFs
        fPlot = fullfile(folder.plots, 'RFs', name, date);
        if ~isfolder(fPlot)
            mkdir(fPlot)
        end
        if isfile(fullfile(f, "_ss_sparseNoise.times.npy"))
            results = io.getNoiseRFData(f);
            for iUnit = 1:length(results.explVars)
                if isnan(results.explVars(iUnit)) || ...
                        results.explVars(iUnit) < minExplainedVariance || ...
                        results.pValues(iUnit) > maxPVal
                    continue
                end
                % rf: [rows x columns x time x ON/OFF];
                rf = squeeze(results.maps(iUnit,:,:,:,:));
                rf(:,:,:,2) = -rf(:,:,:,2);
                [mx,mxTime] = max(max(abs(rf),[],[1 2 4]));
                stimPos = results.edges; % [left, right, top (negative), bottom]
                squW = diff(stimPos(1:2)) / size(rf,1);
                squH = diff(stimPos(3:4)) / size(rf,2);

                figure('Position', [75 195 1470 475])
                for sf = 1:2
                    subplot(1,2,sf)
                    imagesc([stimPos(1)+squW/2 stimPos(2)-squW/2], ...
                        [stimPos(3)+squH/2 stimPos(4)-squH/2], ...
                        rf(:,:,mxTime,sf),[-mx mx])
                    axis image
                    set(gca, 'box', 'off')
                    colormap(gca, cms(:,:,sf))
                    title(titles{sf})
                    colorbar
                end
                sgtitle(sprintf(...
                    'ROI %d (lam: %.0e, t: %.2fs, EV: %.3f, pVal: %.3f)', ...
                    iUnit, results.lambdas(iUnit), ...
                    results.timestamps(mxTime), results.explVars(iUnit), ...
                    results.pValues(iUnit)))

                saveas(gcf, fullfile(fPlots, sprintf('Unit%03d_noise.jpg', iUnit)));
                close gcf
            end
        end
        if isfile(fullfile(f, "_ss_circles.intervals.npy"))
            results = io.getCircleRFData(f);
            for iUnit = 1:length(results.explVars)
                if isnan(results.explVars(iUnit)) || ...
                        results.explVars(iUnit) < minExplainedVariance || ...
                        results.pValues(iUnit) > maxPVal
                    continue
                end
                % rf: [rows x columns x diameters x time x ON/OFF];
                rf = squeeze(results.maps(iUnit,:,:,:,:,:));
                rf(:,:,:,:,2) = -rf(:,:,:,:,2);
                [~,mxD] = max(max(abs(rf),[],[1 2 4 5]));
                [mx,mxTime] = max(max(abs(rf(:,:,mxD,:,:)),[],[1 2 3 5]));
                rf = squeeze(rf(:,:,mxD,mxTime,:));

                figure('Position', [75 195 1470 475])
                for sf = 1:2
                    subplot(1,2,sf)
                    imagesc(results.x([1 end]), results.y([1 end]), ...
                        rf(:,:,sf),[-mx mx])
                    axis image
                    set(gca, 'box', 'off', 'YDir', 'normal')
                    colormap(gca, cms(:,:,sf))
                    title(titles{sf})
                    colorbar
                end
                sgtitle(sprintf(...
                    'ROI %d (diam: %d, lam: %.0e, t: %.2fs, EV: %.3f, pVal: %.3f)', ...
                    iUnit, results.diameters(mxD), ...
                    results.lambdas(iUnit), results.timestamps(mxTime), ...
                    results.explVars(iUnit), results.pValues(iUnit)))

                saveas(gcf, fullfile(fPlots, sprintf('Unit%03d_circle.jpg', iUnit)));
                close gcf
            end
        end
    end
end