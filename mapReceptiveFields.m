function mapReceptiveFields(folder)

%% Parameters
% for correcting baseline drifts of calcium traces at start of experiments
win_decay = 20; % in s, window to test whether baseline is higher than normal
thresh_decay = 1.5; % in std, threshold for drift
win_correct = 150; % in s, window to fit exponential

% to high-pass filter traces
smoothWin = 20; % in s

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
% kernels and RFs simultaneously
lambdas = logspace(-4, -1, 4);
rf_timeLimits = [0.2 0.4];
crossFolds = 10;

% for evaluation of receptive fields (significance/goodness)
minExplainedVariance = 0.01;
maxPVal = 0.05;

% thresholds for receptive fields used for eye shifts
minEV_shift = 0.04;
minPeakNoiseRatio_shift = 7.7;

% fit Gaussian to RF
thresh_diam = 0.75;
thresh_subfield = 0.7;

% for plotting RFs
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
ellipse_x = linspace(-pi, pi, 100);
titles = {'ON field','OFF field'};
RFtypes = {'ON', 'OFF', 'ON+OFF'};

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
        fprintf('  %s %s\n', name, date)
        f = fullfile(folder.data, name, date, '001');

        % ignore datasets for which white noise or circle data do not exist
        if (~isfile(fullfile(f, "_ss_sparseNoise.times.npy")) && ...
                ~isfile(fullfile(f, "_ss_circles.intervals.npy")))
            continue
        end

        folderRes = fullfile(folder.results, name, date);
        if ~isfolder(folderRes)
            mkdir(folderRes);
        end

        %% Load and prepare data
        caData = io.getCalciumData(f);
        for stimType = 1:2
            if stimType == 1
                % if white noise data available, use it to map RFs
                if isfile(fullfile(f, "_ss_sparseNoise.times.npy"))
                    stimData = io.getVisNoiseInfo(f);
                else
                    continue
                end
                % edges: [left right top bottom] (above horizon: >0)
                edges = double([stimData.edges([1 2]), -stimData.edges([3 4])]);
                gridW = diff(edges(1:2)) / size(stimData.frames,3);
                gridH = -diff(edges(3:4)) / size(stimData.frames,2);
                % ignore pixels in ipsilateral (right) hemifield
                if edges(1) * edges(2) < 0
                    % determine right edge of all pixel columns
                    rightEdges = edges(1) + ...
                        (1:size(stimData.frames,3)) .* gridW;
                    validPix = find(rightEdges <= 0);
                    stimData.frames = stimData.frames(:,:,validPix);
                    edges(2) = rightEdges(validPix(end));
                end
                stimMatrix = stimData.frames(stimData.stimOrder,:,:);
                stimSize = size(stimMatrix, [2 3]);
            elseif stimType == 2
                % if circle data available, use that to map RFs
                if isfile(fullfile(f, "_ss_circles.intervals.npy"))
                    stimData = io.getCircleInfo(f);
                else
                    continue
                end
                % create stimulus matrix [t_stim x rows x columns x circle_diameter]
                [stimMatrix, x, y, diameters] = ...
                    circles.getStimMatrix(stimData.times, ...
                    stimData.xyPos, stimData.diameter, stimData.isWhite);
                stimSize = [length(y) length(x) length(diameters)];
                gridW = median(diff(x));
                gridH = median(-diff(y));
                % edges: [left right top bottom] (above horizon: >0)
                edges = [x(1)-0.5*gridW x(end)+0.5*gridW ...
                    y(1)+0.5*gridH y(end)-0.5*gridH];
            end

            t_stim = stimData.times;
            tBin_stim = median(diff(t_stim));
            t_rf = (floor(rf_timeLimits(1) / tBin_stim) : ...
                ceil(rf_timeLimits(2) / tBin_stim)) .* tBin_stim;
            rfBins = t_rf/tBin_stim;

            % generate toeplitz matrix for stimulus
            [toeplitz, t_toeplitz] = ...
                whiteNoise.makeStimToeplitz(stimMatrix, t_stim, rfBins);

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

            % high-pass filter traces: remove smoothed traces
            caTraces = traces.highPassFilter(caTraces, t_ca, smoothWin);

            %% Map RFs
            % rFields: [rows x columns (x diameters) x t_rf x ON/OFF x units]
            [rFields, ev] = ...
                whiteNoise.getReceptiveField(caTraces, t_ca, ...
                toeplitz, t_toeplitz, stimSize, rfBins, ...
                lambdas, crossFolds);

            v = squeeze(mean(ev,3)); % [neuron x lamStim]
            % average EV across cross-folds
            [maxEV, maxStimLam] = max(v,[],2);
            bestLambdas = lambdas(maxStimLam)';

            % test signficance of each RF (note: resulting ev are not
            % cross-validated, while maxEV are)
            [ev, ev_shift] = ...
                whiteNoise.receptiveFieldShiftTest( ...
                caTraces, t_ca, toeplitz, t_toeplitz, ...
                stimSize, rfBins, rFields, bestLambdas, 500);
            pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
            pvals(isnan(ev)) = NaN;

            %--------------------------------------------------------------
            % Uncomment if RFs have already been mapped, and only Gaussian
            % fit should be performed.
            %
            % if stimType == 1
            %     % if white noise data available, use it to map RFs
            %     if isfile(fullfile(f, "_ss_sparseNoise.times.npy"))
            %         results = io.getNoiseRFData(folderRes);
            %     else
            %         continue
            %     end
            %     rFields = results.maps;
            % elseif stimType == 2
            %     % if circle data available, use that to map RFs
            %     if isfile(fullfile(f, "_ss_circles.intervals.npy"))
            %         results = io.getCircleRFData(folderRes);
            %     else
            %         continue
            %     end
            %     rFields = results.maps;
            % end
            % dims = ndims(rFields);
            % rFields = permute(rFields, [2:dims 1]);
            % pvals = results.pValues;
            % maxEV = results.explVars;
            % bestLambdas = results.lambdas;
            % t_rf = results.timestamps;
            %--------------------------------------------------------------

            % fit Gaussian to RF
            % find neurons with good enough RFs to do fit Gaussian
            validRF = find(pvals < maxPVal & maxEV > minExplainedVariance);
            % rfGaussPars: [amplitude, xSTD, yCentre, ySTD, rotation]
            rfGaussPars = NaN(length(maxEV), 6);
            peakNoiseRatio = NaN(length(maxEV), 1);
            bestSubFields = NaN(length(maxEV), 1);
            subFieldSigns = NaN(length(maxEV), 2);
            optimalDelays = NaN(length(maxEV), 1);
            if stimType == 2 % circle paradigm
                sizeTuning = NaN(length(maxEV), length(diameters));
            end
            % perform gaussian fit to determine x-position for valid RFs
            for cellID = validRF'
                if stimType == 1 % visual noise
                    rf = squeeze(rFields(:,:,:,:,cellID));
                    % invert polarity of OFF field so that positive values
                    % mean: unit is driven by black square
                    rf(:,:,:,2) = -rf(:,:,:,2);
                else % circle paradigm
                    % continue with good circle sizes (average across sizes
                    % that drive unit well)
                    rf = squeeze(rFields(:,:,:,:,:,cellID));
                    % invert polarity of OFF field so that positive values
                    % mean: unit is driven by black circle
                    rf(:,:,:,:,2) = -rf(:,:,:,:,2);
                    % average across time
                    rf_tmp = mean(rf,4,"omitnan");
                    % most driven RF "pixel"
                    [~,mpx] = max(rf_tmp, [], "all");
                    % translate "pixel" to indices in RF dimensions
                    [mr,mc,md,~,ms] = ind2sub(size(rf_tmp), mpx);
                    % get size tuning for optimal RF location and subfield
                    mx = squeeze(rf_tmp(mr,mc,:,1,ms));
                    % sign of strongest response (enhanced or suppressed)
                    sgn = sign(mx(md));
                    sizeTuning(cellID,:) = mx;
                    % determine circle sizes that drive cell to 
                    % >thresh_diam of maximum response
                    gd = (mx .* sgn) ./ max(mx.*sgn) > thresh_diam;
                    % average across well-driving circle sizes
                    rf = squeeze(mean(rf(:,:,gd,:,:),3));
                end
                % rf: [rows x columns x delays x subfields]
                % find best subfield (combination)
                % average across time
                rf_tmp = squeeze(mean(rf,3));
                signs = NaN(1,2);
                subs = NaN(1,3);
                % find max response and sign for each subfield
                for sub = 1:2
                    r = rf_tmp(:,:,sub);
                    [subs(sub),ind] = max(abs(r), [], "all");
                    signs(sub) = sign(r(ind));
                end
                % average across subfields (using signs)
                rf_tmp(:,:,3) = (rf_tmp(:,:,1).*signs(1) + ...
                    rf_tmp(:,:,2).*signs(2)) ./ 2;
                subs(3) = max(rf_tmp(:,:,3), [], "all");
                [m, mxSub] = max(subs);
                % use ON and OFF subfields if average is >thresh_subfield *
                % max of best single subfield
                if subs(3) > thresh_subfield * m
                    mxSub = 3;
                end
                bestSubFields(cellID) = mxSub;
                subFieldSigns(cellID,:) = signs;
                % consider only best subfield for fitting
                if mxSub < 3
                    rf = rf(:,:,:,mxSub) .* signs(mxSub);
                else
                    rf = (rf(:,:,:,1) .* signs(1) + ...
                        rf(:,:,:,2) .* signs(2)) ./ 2;
                end
                % find optimal delay for fitting
                % take max across pixels
                [~,mxTime] = max(max(rf, [], [1 2]));
                optimalDelays(cellID) = mxTime;
                rf = rf(:,:,mxTime);

                % interpolate RF so that pixels are square with edge length
                % of 1 visual degree
                % note: xx and yy are locations of rf_visDeg relative to
                % pixel locations of rf (i.e. 1 is at the centre of the 1st
                % pixel of rf)
                [rf_visDeg,xx,yy] = whiteNoise.interpolateRFtoVisualDegrees(...
                    rf, edges);

                % fit 2D Gausian and shift horzontal position of RFs across
                % eye positions
                [fitPars, rf_gauss] = whiteNoise.fit2dGaussRF(rf_visDeg, false);
                % transform RF position to absolute values (relative to
                % visual field)
                fitPars(2) = edges(1) + xx(1)*gridW + fitPars(2)-1;
                fitPars(4) = edges(3) - yy(1)*gridH - (fitPars(4)-1);

                % get peak-to-noise ratio
                noise = std(rf_visDeg - rf_gauss, 0, "all");
                peakNoiseRatio(cellID) = fitPars(1) / noise;

                % collect fitted parameters
                rfGaussPars(cellID,:) = fitPars;
            end

            % save results
            dims = 1:ndims(rFields);
            results.maps = permute(rFields, dims([end 1:end-1]));
            results.explVars = maxEV;
            results.lambdas = bestLambdas;
            results.pValues = pvals;
            results.timestamps = t_rf;
            results.gaussPars = rfGaussPars;
            results.peakToNoise = peakNoiseRatio;
            results.bestSubfields = bestSubFields;
            results.subfieldSigns = subFieldSigns;
            results.optimalDelays = optimalDelays;
            if stimType == 1
                results.edges = edges;
                io.writeNoiseRFResults(results, folderRes);
            elseif stimType == 2
                results.x = x;
                results.y = y;
                results.diameters = diameters;
                results.sizeTuning = sizeTuning;
                io.writeCircleRFResults(results, folderRes);
            end
            clear results
        end

        %% Plot RFs
        fPlots = fullfile(folder.plots, 'RFs', name, date);
        if ~isfolder(fPlots)
            mkdir(fPlots)
            mkdir(fullfile(fPlots, 'doNotShift'))
        end
        if isfile(fullfile(f, "_ss_sparseNoise.times.npy"))
            results = io.getNoiseRFData(folderRes);
            for iUnit = 1:length(results.explVars)
                if isnan(results.explVars(iUnit)) || ...
                        results.explVars(iUnit) <= minExplainedVariance || ...
                        results.pValues(iUnit) >= maxPVal
                    continue
                end
                % rf: [rows x columns x time x ON/OFF];
                rf = squeeze(results.maps(iUnit,:,:,:,:));
                rf(:,:,:,2) = -rf(:,:,:,2);
                mxTime = results.optimalDelays(iUnit);
                mx = max(max(abs(rf(:,:,mxTime,:)),[],"all"));
                edges = results.edges; % [left, right, top, bottom]
                gridW = diff(edges(1:2)) / size(rf,1);
                gridH = -diff(edges(3:4)) / size(rf,2);
                rfGaussPars = results.gaussPars(iUnit,:);

                figure('Position', [75 195 1470 475])
                tiledlayout(2, 1, "TileSpacing", "tight")
                for sf = 1:2
                    nexttile
                    imagesc([edges(1)+gridW/2 edges(2)-gridW/2], ...
                        [edges(3)-gridH/2 edges(4)+gridH/2], ...
                        rf(:,:,mxTime,sf),[-mx mx])
                    hold on
                    % ellipse at 2 STD (x and y), not rotated, not shifted
                    x = rfGaussPars(3) * cos(ellipse_x) * 2;
                    y = rfGaussPars(5) * sin(ellipse_x) * 2;
                    % rotate and shift ellipse
                    x_rot = rfGaussPars(2) + ...
                        x .* cos(rfGaussPars(6)) - ...
                        y .* sin(rfGaussPars(6));
                    y_rot = rfGaussPars(4) + ...
                        x .* sin(rfGaussPars(6)) + ...
                        y .* cos(rfGaussPars(6));
                    n = x_rot < edges(1) | x_rot > edges(2) | ...
                        y_rot > edges(3) | y_rot < edges(4);
                    x_rot(n) = NaN;
                    y_rot(n) = NaN;
                    plot(x_rot, y_rot, 'k')
                    axis image off
                    if sf == 2
                        axis on
                        set(gca, 'box', 'off')
                    end
                    set(gca, 'YDir', 'normal')
                    colormap(gca, cms(:,:,sf))
                    title(titles{sf})
                    colorbar
                end
                sgtitle(sprintf(...
                    'ROI %d: %s (lam: %.0e, t: %.2fs, EV: %.3f, pVal: %.3f, peak/noise: %.1f)', ...
                    iUnit, RFtypes{results.bestSubfields(iUnit)}, ...
                    results.lambdas(iUnit), ...
                    results.timestamps(mxTime), results.explVars(iUnit), ...
                    results.pValues(iUnit), results.peakToNoise(iUnit)))

                if results.explVars(iUnit) > minEV_shift && ...
                        results.peakToNoise(iUnit) > minPeakNoiseRatio_shift
                    saveas(gcf, fullfile(fPlots, ...
                        sprintf('Unit%03d_noise.jpg', iUnit)));
                else
                    saveas(gcf, fullfile(fPlots, 'doNotShift', ...
                        sprintf('Unit%03d_noise.jpg', iUnit)));
                end
                close gcf
            end
        end
        if isfile(fullfile(f, "_ss_circles.intervals.npy"))
            results = io.getCircleRFData(folderRes);
            for iUnit = 1:length(results.explVars)
                if isnan(results.explVars(iUnit)) || ...
                        results.explVars(iUnit) <= minExplainedVariance || ...
                        results.pValues(iUnit) >= maxPVal
                    continue
                end
                mxTime = results.optimalDelays(iUnit);
                % rf: [rows x columns x diameters x time x ON/OFF];
                rf = squeeze(results.maps(iUnit,:,:,:,mxTime,:));
                rf(:,:,:,2) = -rf(:,:,:,2);
                % determine best delay; average across diameters
                % r = mean(rf,3,"omitnan");
                % [~,mxTime] = max(max(abs(r),[],[1 2 3 5]));
                % rf = squeeze(rf(:,:,:,mxTime,:));
                mx = max(abs(rf),[],"all");
                diams = results.diameters;
                gridW = median(diff(results.x));
                gridH = median(-diff(results.y));
                % edges: [left right top bottom] (above horizon: >0)
                edges = [results.x(1)-0.5*gridW results.x(end)+0.5*gridW ...
                    results.y(1)+0.5*gridH results.y(end)-0.5*gridH];
                rfGaussPars = results.gaussPars(iUnit,:);

                figure('Position', [75 195 1470 475])
                tiledlayout(2, length(diams), "TileSpacing", "tight");
                for sf = 1:2
                    for d = 1:length(diams)
                        nexttile
                        imagesc(results.x([1 end]), results.y([1 end]), ...
                            rf(:,:,d,sf),[-mx mx])
                        hold on
                        % ellipse at 2 STD (x and y), not rotated, not shifted
                        x = rfGaussPars(3) * cos(ellipse_x) * 2;
                        y = rfGaussPars(5) * sin(ellipse_x) * 2;
                        % rotate and shift ellipse
                        x_rot = rfGaussPars(2) + ...
                            x .* cos(rfGaussPars(6)) - ...
                            y .* sin(rfGaussPars(6));
                        y_rot = rfGaussPars(4) + ...
                            x .* sin(rfGaussPars(6)) + ...
                            y .* cos(rfGaussPars(6));
                        n = x_rot < edges(1) | x_rot > edges(2) | ...
                            y_rot > edges(3) | y_rot < edges(4);
                        x_rot(n) = NaN;
                        y_rot(n) = NaN;
                        plot(x_rot, y_rot, 'k')
                        axis image off
                        set(gca, 'box', 'off', 'YDir', 'normal')
                        colormap(gca, cms(:,:,sf))
                        if sf==1
                            title(sprintf('%.1f deg diam.', diams(d)))
                        end
                        if d==length(diams)
                            axis on
                            set(gca, 'box', 'off')
                            colorbar
                        end
                    end
                end
                gd = results.sizeTuning(iUnit,:);
                [~,md] = max(abs(gd));
                sgn = sign(gd(md));
                gd = (gd .* sgn) ./ max(gd.*sgn) > thresh_diam;
                sgtitle(sprintf(...
                    'ROI %d: %s (lam: %.0e, t: %.2fs, EV: %.3f, pVal: %.3f, peak/noise: %.1f, diameters:%s)', ...
                    iUnit, RFtypes{results.bestSubfields(iUnit)}, ...
                    results.lambdas(iUnit), ...
                    results.timestamps(mxTime), results.explVars(iUnit), ...
                    results.pValues(iUnit), results.peakToNoise(iUnit), ...
                    sprintf(' %.1f',diams(gd))))

                if results.explVars(iUnit) > minEV_shift && ...
                        results.peakToNoise(iUnit) > minPeakNoiseRatio_shift
                    saveas(gcf, fullfile(fPlots, ...
                        sprintf('Unit%03d_circle.jpg', iUnit)));
                else
                    saveas(gcf, fullfile(fPlots, 'doNotShift', ...
                        sprintf('Unit%03d_circle.jpg', iUnit)));
                end
                close gcf
            end
        end
    end
end