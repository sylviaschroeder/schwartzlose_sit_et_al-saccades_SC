%% Parameters
% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = 0.01;
minPVal = 0.05;
maxLambda = 1;
minEVnew = 0.2;

% for binning of horizontal eye position
numEyeBins = 10;

% for correcting baseline drifts of calcium traces at start of experiments
driftWin = 20; % in s, window to test whether baseline is higher than normal
driftThresh = 1.5; % in std, threshold for drift
correctWin = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously)
lambdasStim = 0.001;
RFlimits = [0.2 0.4];
crossFolds = 1;
RFtypes = {'ON', 'OFF', 'ON+OFF'};

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
        t_stim = noiseData.times;
        eyePos = interp1(pupilData.time, pupilData.pos(:,1), t_stim);

        % percentiles of horizontal eye positions (including data from
        % the whole session)
        % Note: position from nasal (smaller) to temporal (larger)
        eyePosQuart = prctile(pupilData.pos(:,1), [25 50 75]);
        eyePosNorm = (eyePos - eyePosQuart(2)) ./ diff(eyePosQuart([1 3]));
        [~,eyeBinEdges,bin] = histcounts(eyePosNorm, ...
            prctile(eyePosNorm, linspace(0,100,numEyeBins+1)));

        % find neurons with good RF, ignore the rest
        validRF = find(rfData.pValues < minPVal & ...
            rfData.explVars > minExplainedVarianceStim & ...
            rfData.lambdas < maxLambda);

        % interpolate calcium traces to align all to same time
        t_ind = caData.time > t_stim(1) - 10 & caData.time < t_stim(end) + 10;
        traces = caData.traces(t_ind,validRF);
        t_ca = caData.time(t_ind);
        timeBin = median(diff(t_ca));
        if length(caData.delays) > 1
            delta_t = median(diff(caData.delays));
            upsample = round(timeBin / delta_t);
            timeBin = timeBin / upsample;
            t_up = reshape((t_ca + (0:upsample-1) * timeBin)', [], 1);
            tr_up = NaN(length(t_up), size(traces,2));
            for d = 1:length(caData.delays)
                indUnits = find(caData.planes(validRF) == d);
                for n = indUnits'
                    if all(isnan(traces(:,n)))
                        continue
                    end
                    nanInd1 = isnan(traces(:,n));
                    tr_up(:,n) = interp1(t_ca(~nanInd1) + caData.delays(d), ...
                        traces(~nanInd1,n), t_up, 'pchip');
                    nanInd2 = reshape(repmat(nanInd1, 1, upsample), [], 1);
                    tr_up(nanInd2,n) = NaN;
                end
            end
            t_ca = t_up;
            traces = tr_up;
        end
        
        % remove strong baseline decay at start of experiment in cells that
        % show it
        indUnits = find(mean(traces(1:round(driftWin / timeBin),:), 1, 'omitnan') > ...
            mean(traces, 1, 'omitnan') + driftThresh .* std(traces,0,1, 'omitnan'));
        ind = round(correctWin / timeBin);
        for iUnit = 1:length(indUnits)
            y = traces(:, indUnits(iUnit));
            y = fillmissing(y, 'linear');
            % fit double exponential to start of trace
            f = fit((1:length(y))', y, ...
                @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
                'Lower', [0 0 0 0 0], ...
                'Upper', [max(y) max(y) 500 max(y) 500], ...
                'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
            % remove fit
            traces(:, indUnits(iUnit)) = y - f(1 : size(traces,1)) + f.a;
        end

        % map RFs for each eye position
        frames = noiseData.frames(noiseData.stimOrder,:,:);
        stimFrameDur = median(diff(t_stim));
        RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
            ceil(RFlimits(2) / stimFrameDur);
        RFs = NaN(size(frames,2), size(frames,3), length(RFtimesInFrames), ...
            2, length(validRF), numEyeBins);
        ev = NaN(length(validRF), numEyeBins);
        for p = 1:numEyeBins
            fr = frames;
            % set all stimulus frames where eye was outside the position of
            % interest to NaN so they don't contribute to the RF fitting
            fr(bin ~= p,:,:) = NaN;
            [RFs(:,:,:,:,:,p), ev(:,p)] = ...
                whiteNoise.getReceptiveField(traces, t_ca, fr, ...
                t_stim, RFtimesInFrames, lambdasStim, crossFolds);
        end
        % positive pixel in OFF subfield means: driven by black
        RFs(:,:,:,2,:,:) = -RFs(:,:,:,2,:,:);
        ev = mean(ev,2);

        % choose neurons with highest explained variances
        goodUnits = find(ev > minEVnew);

        for iUnit = 1:length(goodUnits)
            % find best subfield (ON or OFF or ON/OFF), find best time of RF
            rf = mean(RFs(:,:,:,:,goodUnits(iUnit),:),6);
            [~,mxTime] = max(max(reshape(permute(abs(rf), [1 2 4 3]), ...
                [], size(RFs,3)), [], 1));
            rf = squeeze(rf(:,:,mxTime,:));
            signs = NaN(1,2);
            subs = NaN(1,3);
            for sub = 1:2
                r = rf(:,:,sub);
                [subs(sub),ind] = max(abs(r), [], "all");
                signs(sub) = sign(r(ind));
            end
            subs(3) = max(mean(abs(rf),3), [], "all");
            [~,mxSub] = max(subs);
            if mxSub < 3
                rf = squeeze(RFs(:,:,mxTime,mxSub,goodUnits(iUnit),:)) .* signs(mxSub);
            else
                rf = squeeze(RFs(:,:,mxTime,1,goodUnits(iUnit),:) .* signs(1) + ...
                    RFs(:,:,mxTime,2,goodUnits(iUnit),:) .* signs(2)) ./ 2;
            end

            % interpolate RF so that pixels are square with edge length of 
            % 1 visual degree
            rf_visDeg = whiteNoise.interpolateRFtoVisualDegrees(rf, noiseData.edges);

            % fit 2D Gausian and shift parameter to RFs across eye positions
            [fitPars, rf_fitted] = whiteNoise.fit2dGaussRF_wShift(rf_visDeg, false);
        end
    end
end