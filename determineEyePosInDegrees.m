function determineEyePosInDegrees(folder)

%% Parameters
% thresholds for all receptive fields
minEV_shift = 0.01; % increase to 0.05 for shifts?
maxPVal = 0.05;
% thresholds for receptive fields used for eye shifts
minPeakNoiseRatio_shift = 7.7;
% thresholds for rest of receptive fields (only determine median position)
minPeakNoiseRatio_median = 5;

% for binning of horizontal eye position
numEyeBinsNoise = 9;
numEyeBinsCircles = 5;

% for correcting baseline drifts of calcium traces at start of experiments
driftWin = 20; % in s, window to test whether baseline is higher than normal
driftThresh = 1.5; % in std, threshold for drift
correctWin = 150; % in s, window to fit exponential

% to high-pass filter traces
smoothWin = 20; % in s

% for receptive field mapping
thresh_diam = 0.75;
lambdasStim = 0.002; %[0.002 0.1];
rf_timeLimits = [0.2 0.4];
crossFolds = 1;
RFtypes = {'ON', 'OFF', 'ON+OFF'};
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
ellipse_x = linspace(-pi, pi, 100);

%% Loop through all datasets
% Note: for all SC neurons recorded with 2P, recording was in right
% hemisphere and left eye was on video
subjects = dir(fullfile(folder.data));
subjects = subjects(~startsWith({subjects.name},'.'));

for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        f = fullfile(folder.data, name, date, '001');
        folderRes = fullfile(folder.results, name, date);

        if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy')) && ...
                ~isfile(fullfile(f, "_ss_circles.intervals.npy"))
            continue
        end

        fprintf('  %s %s\n', name, date)
        folderPl = fullfile(folder.plots, 'RFsByEyePos', name, date);
        if ~isfolder(folderPl)
            mkdir(folderPl)
        end
        
        %% Load and prepare data
        caData = io.getCalciumData(f);
        pupilData = io.getPupilData(f);
        if isfile(fullfile(folderRes, "_ss_circlesRf.pValues.npy"))
            stimData = io.getCircleInfo(f);
            rfData = io.getCircleRFData(folderRes);
            
            gridW = median(diff(rfData.x));
            gridH = median(-diff(rfData.y));
            % edges: [left right top bottom] (above horizon: >0)
            edges = [rfData.x(1)-0.5*gridW rfData.x(end)+0.5*gridW ...
                rfData.y(1)+0.5*gridH rfData.y(end)-0.5*gridH];

            nb = numEyeBinsCircles;
        elseif isfile(fullfile(f, "_ss_rf.pValues.npy"))
            stimData = io.getVisNoiseInfo(f);
            rfData = io.getNoiseRFData(folderRes);

            % edges: [left right top bottom] (above horizon: >0)
            edges = rfData.edges;
            gridW = diff(edges(1:2)) / size(stimData.frames,3);
            gridH = -diff(edges(3:4)) / size(stimData.frames,2);
            validPix = 1:size(rfData.maps, 3);
            if stimData.edges(1) * stimData.edges(2) < 0
                stimData.frames = stimData.frames(:,:,validPix);
            end

            nb = numEyeBinsNoise;
        else
            continue
        end
        
        t_stim = stimData.times;
        tBin_stim = median(diff(t_stim));
        t_rf = (floor(rf_timeLimits(1) / tBin_stim) : ...
            ceil(rf_timeLimits(2) / tBin_stim)) .* tBin_stim;
        rfBins = t_rf/tBin_stim;

        if isfield(stimData, "diameter") % circle paradigm
            % create stimulus matrix [t_stim x rows x columns x circle_diameter]
            [stimMatrix, x, y, diameters] = circles.getStimMatrix(t_stim, ...
                stimData.xyPos, stimData.diameter, stimData.isWhite);
            stimSize = [length(y) length(x) length(diameters)];
        else % visual noise
            stimMatrix = stimData.frames(stimData.stimOrder,:,:);
            stimSize = size(stimMatrix, [2 3]);
        end
        % generate toeplitz matrix for stimulus
        [toeplitz, t_toeplitz] = ...
            whiteNoise.makeStimToeplitz(stimMatrix, t_stim, rfBins);

        %% Prepare calcium traces
        validRF = find(rfData.explVars > minEV_shift & ...
            rfData.pValues < maxPVal);
        % interpolate calcium traces to align all to same time
        t_ind = caData.time > t_stim(1) - 10 & caData.time < t_stim(end) + 10;
        caTraces = caData.traces(t_ind,validRF);
        t_ca = caData.time(t_ind);
        [caTraces, t_ca] = traces.alignSampling(caTraces, t_ca, ...
            caData.planes(validRF), caData.delays);

        % remove strong baseline decay at start of experiment in cells that
        % show it
        caTraces = traces.removeDecay(caTraces, t_ca, driftWin, ...
            correctWin, driftThresh);

        % high-pass filter traces: remove smoothed traces
        caTraces = traces.highPassFilter(caTraces, t_ca, smoothWin);

        %% Map shifted RFs
        % find neurons with good enough RFs to do shifts
        validForShift = find(rfData.explVars > minEV_shift & ...
            rfData.pValues < maxPVal & ...
            rfData.peakToNoise > minPeakNoiseRatio_shift);
        indValidShift = ismember(validRF, validForShift);

        % interpolate pupil data to match stimulus times
        tBin_eye = median(diff(pupilData.time), "omitnan");
        down = ceil(tBin_stim / tBin_eye / 2) * 2 + 1;
        eyePos = medfilt1(pupilData.pos(:,1), down, "omitnan");
        nanInd1 = isnan(eyePos);
        eyePos = interp1(pupilData.time(~nanInd1), ...
            eyePos(~nanInd1), t_toeplitz);
        nanInd2 = histcounts(pupilData.time(nanInd1), t_toeplitz) > 0;
        eyePos(nanInd2) = NaN;
        [~,~,bin] = histcounts(eyePos, ...
            prctile(eyePos, linspace(0, 100, nb+1)));
        % get median eye position for each bin
        eyePosPerBin = NaN(1, nb);
        for b = 1:nb
            eyePosPerBin(b) = median(eyePos(bin == b), "omitnan");
        end

        % RFs_shift: [units x rows x columns (x diameters) x delays x
        % subfields x shifts]
        RFs_shift = NaN([length(validForShift) stimSize length(t_rf) 2 0]);
        numDim = length(stimSize) + 3;
        for p = 1:nb
            % ignore all time points for which eye position was not in bin p
            ignore = bin ~= p;
            rfs = whiteNoise.getReceptiveField(caTraces(:,indValidShift), t_ca, ...
                toeplitz, t_toeplitz, stimSize, rfBins, ...
                lambdasStim, crossFolds, ignore);
            RFs_shift = cat(numDim+1, RFs_shift, permute(rfs, [numDim 1:numDim-1]));
        end
        if isfield(stimData, "diameter")
            % continue with average across good circle sizes
            for iUnit = 1:size(RFs_shift,1)
                cellID = validForShift(iUnit);
                gd = rfData.sizeTuning(cellID,:);
                [~,md] = max(abs(gd));
                sgn = sign(gd(md));
                gd = (gd .* sgn) ./ max(gd.*sgn) > thresh_diam;
                RFs_shift(iUnit,:,:,1,:,:,:) = ...
                    mean(RFs_shift(iUnit,:,:,gd,:,:,:),4,"omitnan");
            end
            RFs_shift = permute(RFs_shift(:,:,:,1,:,:,:),[1 2 3 5 6 7 4]);
        end
        % invert polarity of OFF field so that positive values mean: unit 
        % is driven by black square/circle
        RFs_shift(:,:,:,:,2,:) = -RFs_shift(:,:,:,:,2,:);

        %% Map median position RFs
        % find neurons where RF was not good enough to do shifts
        validNoShift = validRF(~indValidShift);
        % map RFs for eye positions within quartiles
        [~,~,bin] = histcounts(eyePos, prctile(eyePos, [25 75]));

        % ignore all time points for which eye position was within quartiles
        % RFs_median: [rows x columns (x diameters) x delays x subfields x units]
        RFs_median = whiteNoise.getReceptiveField( ...
            caTraces(:,~indValidShift), t_ca, ...
            toeplitz, t_toeplitz, stimSize, rfBins, ...
            lambdasStim, crossFolds, bin<1);
        numDim = length(stimSize) + 3;
        % RFs_median: [units x rows x columns (x diameters) x delays x subfields]
        RFs_median = permute(RFs_median, [numDim 1:numDim-1]);
        if isfield(stimData, "diameter")
            % continue with average across good circle sizes
            for iUnit = 1:size(RFs_median,1)
                cellID = validNoShift(iUnit);
                gd = rfData.sizeTuning(cellID,:);
                [~,md] = max(abs(gd));
                sgn = sign(gd(md));
                gd = (gd .* sgn) ./ max(gd.*sgn) > thresh_diam;
                RFs_median(iUnit,:,:,1,:,:) = ...
                    mean(RFs_median(iUnit,:,:,gd,:,:),4,"omitnan");
            end
            RFs_median = permute(RFs_median(:,:,:,1,:,:),[1 2 3 5 6 4]);
        end
        % invert polarity of OFF field so that positive values mean: unit 
        % is driven by black square/circle
        RFs_median(:,:,:,:,2) = -RFs_median(:,:,:,:,2);

        %% Fit shifted RF maps with 2D Gaussian
        % rfGaussShift: [units x eye position], horizontal position of RF
        % for each eye position
        rfGaussShift = NaN(length(rfData.explVars), nb);
        % rfGaussConst: [units x Gauss parameters (amplitude, xSTD, 
        % yCentre, ySTD, rotation)]
        rfGaussConst = NaN(length(rfData.explVars), 5);
        for iUnit = 1:length(validForShift)
            cellID = validForShift(iUnit);
            mxSub = rfData.bestSubfields(cellID);
            signs = rfData.subfieldSigns(cellID,:);
            mxTime = rfData.optimalDelays(cellID);

            % only consider best subfield (or average)
            if mxSub < 3
                rf = squeeze(RFs_shift(iUnit,:,:,:,mxSub,:)) .* signs(mxSub);
            else
                rf = squeeze(RFs_shift(iUnit,:,:,:,:,:));
                rf(:,:,:,1,:) = rf(:,:,:,1,:) .* signs(1);
                rf(:,:,:,2,:) = rf(:,:,:,2,:) .* signs(2);
                rf = squeeze(mean(rf, 4, "omitnan"));
            end
            % only consider best delay
            rf = squeeze(rf(:,:,mxTime,:));
            rf(isnan(rf)) = 0;

            % interpolate RF so that pixels are square with edge length of 
            % 1 visual degree
            % note: xx and yy are locations of rf_visDeg relative to pixel
            % locations of rf (i.e. 1 is at the centre of the 1st pixel of
            % rf)
            [rf_visDeg,xx,yy] = whiteNoise.interpolateRFtoVisualDegrees(...
                rf, edges);

            % fit 2D Gausian and shift parameter to RFs across eye positions
            fitPars = whiteNoise.fit2dGaussRF_wShift(rf_visDeg, false);
            % transform RF position to absolute values (relative to visual
            % field)
            fitPars(2) = edges(1) + xx(1)*gridW + fitPars(2)-1;
            fitPars(4) = edges(3) - yy(1)*gridH - (fitPars(4)-1);

            % collect fitted parameters
            rfGaussConst(cellID,:) = fitPars([1 3:6]);
            xCentres = fitPars(2) + [0 fitPars(7:end)];
            rfGaussShift(cellID,:) = xCentres;
        end

        %% Fit median RF maps with 2D Gaussian
        medianRFpositions = NaN(length(rfData.explVars), 1);
        peakNoiseRatio = NaN(length(rfData.explVars), 1);
        for iUnit = 1:length(validNoShift)
            cellID = validNoShift(iUnit);
            mxSub = rfData.bestSubfields(cellID);
            signs = rfData.subfieldSigns(cellID,:);
            mxTime = rfData.optimalDelays(cellID);

             % only consider best subfield (or average)
            if mxSub < 3
                rf = squeeze(RFs_median(iUnit,:,:,:,mxSub)) .* signs(mxSub);
            else
                rf = squeeze(RFs_median(iUnit,:,:,:,:));
                rf(:,:,:,1,:) = rf(:,:,:,1) .* signs(1);
                rf(:,:,:,2,:) = rf(:,:,:,2) .* signs(2);
                rf = squeeze(mean(rf, 4, "omitnan"));
            end
            % only consider best delay
            rf = squeeze(rf(:,:,mxTime));
            rf(isnan(rf)) = 0;

            % interpolate RF so that pixels are square with edge length of 
            % 1 visual degree
            % note: xx and yy are locations of rf_visDeg relative to pixel
            % locations of rf (i.e. 1 is at the centre of the 1st pixel of
            % rf)
            [rf_visDeg,xx,yy] = whiteNoise.interpolateRFtoVisualDegrees(...
                rf, edges);

            % fit 2D Gausian and shift parameter to RFs across eye positions
            [fitPars, rf_gauss] = whiteNoise.fit2dGaussRF(rf_visDeg, false);
            % transform RF position to absolute values (relative to visual
            % field)
            fitPars(2) = edges(1) + xx(1)*gridW + fitPars(2)-1;
            fitPars(4) = edges(3) - yy(1)*gridH - (fitPars(4)-1);

            % get peak-to-noise ratio
            noise = std(rf_visDeg - rf_gauss, 0, "all");
            peakNoiseRatio(cellID) = fitPars(1) / noise;

            % collect fitted parameters
            rfGaussConst(cellID,:) = fitPars([1 3:6]);
            medianRFpositions(cellID,:) = fitPars(2);
        end

        %% Linear regression
        % linear regression relating eye position in pixels to RF centre
        % in visual degrees
        shifts = rfGaussShift; % [units x eyePosBins]
        % disregard RF positions that are one RF STD away from stimulus
        % edges
        minSTD = min(rfGaussConst(:,[2 4]),[],2);
        shifts(shifts+minSTD < edges(1) | shifts-minSTD > edges(2)) = NaN;
        % determine median position of RF by taking median across 3 centre
        % positions of eye
        centrePos = [-1 0 1] + ceil(nb/2);
        centrePosValid = find(all(~isnan(shifts(:,centrePos)), 2));
        medianRFpositions(centrePosValid) = ...
            median(shifts(centrePosValid,centrePos), 2);

        if ~isempty(centrePosValid)
            % fit: rf_xPos_norm = a + b * eye_xPos
            % normalize shifted RF positions by each unit's median RF
            % position
            posNorm = shifts(centrePosValid,:) - ...
                medianRFpositions(centrePosValid);
            mdl = fitlm(reshape(repmat(eyePosPerBin, ...
                length(centrePosValid), 1), [], 1), ...
                reshape(posNorm, [], 1), "RobustOpts", "on");
            coeffs = mdl.Coefficients.Estimate;

            % if some shifted RF positions were outside stimulus, use
            % linear model to estimate median RF position based on shifted RF
            % positions inside stimulus:
            % rf_xPos - rf_medianXPos = a + b * eye_xPos
            % => rf_medianXPos = rf_xPos - (a + b * eye_xPos)
            someShiftsValid = setdiff(validForShift, centrePosValid);
            if ~isempty(someShiftsValid)
                invalid = [];
                for cellID = someShiftsValid'
                    x = shifts(cellID,:);
                    ind = ~isnan(x);
                    if all(ind(centrePos))
                        indMed = centrePos;
                    elseif ind(centrePos(2))
                        indMed = centrePos(2);
                    elseif all(ind(centrePos([1 3])))
                        indMed = centrePos([1 3]);
                    else
                        invalid = [invalid; cellID];
                        continue
                    end
                    medianRFpositions(cellID) = ...
                        median(x(indMed)' - predict(mdl, eyePosPerBin(indMed)'));
                end
                someShiftsValid = setdiff(someShiftsValid, invalid);
            end
        else
            coeffs = [NaN NaN];
            someShiftsValid = [];
        end

        %% Save results

        % save RF fits
        rfMapsShifted = NaN([size(rfData.maps,1), size(RFs_shift,2:6)]);
        rfMapsShifted(validForShift,:,:,:,:,:) = RFs_shift;
        writeNPY(rfMapsShifted, fullfile(folderRes, 'rfPerEyePos.maps.npy'))
        rfMaps = NaN([size(rfData.maps,1), size(RFs_median,2:5)]);
        rfMaps(validNoShift,:,:,:,:,:) = RFs_median;
        writeNPY(rfMaps, fullfile(folderRes, 'rfPerEyePos.maps_medianEyePos.npy'))
        writeNPY(rfGaussConst, fullfile(folderRes, 'rfPerEyePos.gaussFitPars_fixed.npy'))
        writeNPY(rfGaussShift, fullfile(folderRes, 'rfPerEyePos.gaussFitPars_xShifts.npy'))
        writeNPY(peakNoiseRatio, fullfile(folderRes, 'rfPerEyePos.peakToNoiseRatio_medianEyePos.npy'))

        % save median horizontal RF positions
        writeNPY(medianRFpositions, fullfile(folderRes, 'rfPerEyePos.medianHorizRFPosition.npy'))
        isMeasured = false(length(caData.ids),1);
        isMeasured(centrePosValid) = true;
        writeNPY(isMeasured, fullfile(folderRes, 'rfPerEyePos.isRFMedian_shift.npy'))
        isModelled = false(length(caData.ids),1);
        isModelled(someShiftsValid) = true;
        writeNPY(isModelled, fullfile(folderRes, 'rfPerEyePos.isRFMedian_modelled.npy'))
        isAcrossAllPupilPos = false(length(caData.ids),1);
        isAcrossAllPupilPos(validNoShift) = true;
        writeNPY(isAcrossAllPupilPos, fullfile(folderRes, 'rfPerEyePos.isRFMedian_medianEyePos.npy'))

        % save median eye positions
        writeNPY(eyePosPerBin, fullfile(folderRes, 'eyeToRFPos.eyeHorizPositions.npy'))

        % save linear regression fit
        writeNPY(coeffs, fullfile(folderRes, 'eyeToRFPos.regressionCoeffs.npy'))
        
        %% Make plots
        % plot RF maps and fitted ellipse for each neuron with shift
        for cellID = validRF'
            mxSub = rfData.bestSubfields(cellID);
            mxTime = rfData.optimalDelays(cellID);
            if isfield(stimData, "diameter")
                gd = rfData.sizeTuning(cellID,:);
                [~,md] = max(abs(gd));
                sgn = sign(gd(md));
                gd = (gd .* sgn) ./ max(gd.*sgn) > thresh_diam;
            end

            if ismember(cellID, validForShift)
                ind = validForShift == cellID;
                rf = squeeze(RFs_shift(ind,:,:,mxTime,:,:));
                mx = max(abs(rf(:)));

                figure('Position', [4 300 1460 420])
                tiledlayout(2, nb, "TileSpacing", "tight")
                for sub = 1:2
                    for pos = 1:nb
                        nexttile
                        imagesc(...
                            [edges(1)+gridW/2 edges(2)-gridW/2], ...
                            [edges(3)-gridH/2 edges(4)+gridH/2], ...
                            rf(:,:,sub,pos),[-mx mx])
                        hold on
                        % ellipse at 2 STD (x and y), not rotated, not shifted
                        x = rfGaussConst(cellID, 2) * cos(ellipse_x) * 2;
                        y = rfGaussConst(cellID, 4) * sin(ellipse_x) * 2;
                        % rotate and shift ellipse
                        x_rot = rfGaussShift(cellID, pos) + ...
                            x .* cos(rfGaussConst(cellID, 5)) - ...
                            y .* sin(rfGaussConst(cellID, 5));
                        y_rot = rfGaussConst(cellID, 3) + ...
                            x .* sin(rfGaussConst(cellID, 5)) + ...
                            y .* cos(rfGaussConst(cellID, 5));
                        n = x_rot < edges(1) | x_rot > edges(2) | ...
                            y_rot > edges(3) | y_rot < edges(4);
                        x_rot(n) = NaN;
                        y_rot(n) = NaN;
                        plot(x_rot, y_rot, 'k')
                        axis image off
                        if sub == 1
                            if pos == 1
                                title('nasal')
                            elseif pos == nb
                                title('temporal')
                            end
                        elseif sub == 2 && pos == nb
                            axis on
                            set(gca, 'box', 'off')
                        end
                        set(gca, 'YDir', 'normal')
                        colormap(gca, cms(:,:,sub))
                    end
                    colorbar
                end
                if isfield(stimData, "diameter")
                    sgtitle(sprintf('Neuron %d: %s (w/o shift: EV: %.3f, peak/noise: %.1f, diameters:%s)', ...
                        cellID, ...
                        RFtypes{mxSub}, rfData.explVars(cellID), ...
                        rfData.peakToNoise(cellID), ...
                        sprintf(' %.1f',diameters(gd))),'FontWeight', 'bold')
                else
                    sgtitle(sprintf('Neuron %d: %s (w/o shift: EV: %.3f, peak/noise: %.1f)', ...
                        cellID, ...
                        RFtypes{mxSub}, rfData.explVars(cellID), ...
                        rfData.peakToNoise(cellID)), 'FontWeight', 'bold')
                end
            else
                ind = validNoShift == cellID;
                rf = squeeze(RFs_median(ind,:,:,mxTime,:));
                mx = max(abs(rf(:)));

                figure('Position', [400 300 700 420])
                tiledlayout(2, 1, "TileSpacing", "tight")
                for sub = 1:2
                    nexttile
                    imagesc(...
                        [edges(1)+gridW/2 edges(2)-gridW/2], ...
                        [edges(3)-gridH/2 edges(4)+gridH/2], ...
                        rf(:,:,sub),[-mx mx])
                    hold on
                    % ellipse at 2 STD (x and y), not rotated, not shifted
                    x = rfGaussConst(cellID, 2) * cos(ellipse_x) * 2;
                    y = rfGaussConst(cellID, 4) * sin(ellipse_x) * 2;
                    % rotate and shift ellipse
                    x_rot = medianRFpositions(cellID) + ...
                        x .* cos(rfGaussConst(cellID, 5)) - ...
                        y .* sin(rfGaussConst(cellID, 5));
                    y_rot = rfGaussConst(cellID, 3) + ...
                        x .* sin(rfGaussConst(cellID, 5)) + ...
                        y .* cos(rfGaussConst(cellID, 5));
                    n = x_rot < edges(1) | x_rot > edges(2) | ...
                        y_rot > edges(3) | y_rot < edges(4);
                    x_rot(n) = NaN;
                    y_rot(n) = NaN;
                    plot(x_rot, y_rot, 'k')
                    axis image off
                    if sub == 2
                        axis on
                        set(gca, 'box', 'off')
                    end
                    set(gca, 'YDir', 'normal')
                    colormap(gca, cms(:,:,sub))
                    colorbar
                end
                if isfield(stimData, "diameter")
                    sgtitle(sprintf('Neuron %d: %s (EV (all eye pos): %.3f, peak/noise: %.1f, diameters:%s)', ...
                        cellID, ...
                        RFtypes{mxSub}, rfData.explVars(cellID), ...
                        peakNoiseRatio(cellID), ...
                        sprintf(' %.1f',diameters(gd))),'FontWeight', 'bold')
                else
                    sgtitle(sprintf('Neuron %d: %s (EV (all eye pos): %.3f, peak/noise: %.1f)', ...
                        cellID, ...
                        RFtypes{mxSub}, rfData.explVars(cellID), ...
                        peakNoiseRatio(cellID)), 'FontWeight', 'bold')
                end
            end
            saveas(gcf, fullfile(folderPl, sprintf('Neuron%03d.jpg', cellID)))
            close gcf
        end
        
        % plot eye position versus horizontal RF centre
        % for mean across neurons: 
        % exclude any neurons where horizontal shift of RF falls outside
        % stimulus edges
        figure('Position', [680 460 990 420])
        tiledlayout(1, 2)
        nexttile
        plot(eyePosPerBin, shifts, 'Color', [1 1 1].*0.5)
        set(gca, 'box', 'off')
        axis tight
        xlabel('Pupil position (pixels)')
        ylabel('Horiz. RF centre (vis deg)')
        title('Absolute RF positions')

        if ~isempty(centrePosValid)
            nexttile
            s = [posNorm; shifts(someShiftsValid,:)-medianRFpositions(someShiftsValid)];
            plot(eyePosPerBin, s', 'Color', [1 1 1] .* 0.5)
            hold on
            plot(eyePosPerBin, median(s, 1, 'omitnan'), 'k', 'LineWidth', 2)
            x = eyePosPerBin([1 end]);
            plot(x, coeffs(1) + x .* coeffs(2), 'r', 'LineWidth', 2)
            set(gca, 'box', 'off')
            axis tight
            xlabel('Pupil position (pixels)')
            ylabel('Horiz. RF centre (median subtracted)')
            title(sprintf('y = %.2f * x + %.2f', coeffs(2), coeffs(1)))
        end
        sgtitle(sprintf('%s %s', name, date), 'FontWeight', 'bold')
        saveas(gcf, fullfile(folder.plots, 'RFsByEyePos', sprintf('%s_%s.jpg', ...
            name, date)))
        close gcf
    end
end
