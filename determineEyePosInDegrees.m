function determineEyePosInDegrees(folder)
%% Parameters
% for evaluation of receptive fields (significance/goodness)
minExplainedVarianceStim = [0.016 0.01];
minPVal = [0.05 0.002];
lambdas = [0.01 0.1];
minEVOld = 0.01;
minPOld = 0.05;
lambdaOld = 1;

% for binning of horizontal eye position
numEyeBins = [10 5];

% for correcting baseline drifts of calcium traces at start of experiments
driftWin = 20; % in s, window to test whether baseline is higher than normal
driftThresh = 1.5; % in std, threshold for drift
correctWin = 150; % in s, window to fit exponential

% for receptive field estimates
% used for fitting 2 RFs (ON and OFF simultaneously)
lambdasStim = [0.002 0.1];
RFlimits = [0.2 0.4];
crossFolds = 1;
RFtypes = {'ON', 'OFF', 'ON+OFF'};
[cm_ON, cm_OFF] = colmaps.getRFMaps;
cms = cat(3, cm_ON, cm_OFF);
ellipse_x = linspace(-pi, pi, 100);

% datasets with lower quality of RFs -> use 2nd set of parameters
lowQualData = {'SS047', {'2015-11-23', '2015-12-03'}; ...
           'SS048', {'2015-12-02'}};

%% Loop through all datasets
% Note: for all SC neurons recorded with 2P, recording was in right
% hemisphere and left eye was on video
fprintf('Fit RFs across pupil positions:\nDatasets:\n')
subjects = dir(fullfile(folder.data, 'SS*'));
for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        f = fullfile(folder.data, name, date, '001');

        if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
            continue
        end
        
        paramSet = 1;
        nb = numEyeBins(1);
        ind = find(strcmp(name, lowQualData(:,1)));
        if ~isempty(ind) && any(strcmp(date, lowQualData{ind,2}))
            paramSet = 2;
            nb = numEyeBins(2);
        end

        fprintf('  %s %s\n', name, date)
        folderPl = fullfile(folder.plots, 'RFsByEyePos', name, date);
        if ~isfolder(folderPl)
            mkdir(folderPl)
        end
        
        % load data
        caData = io.getCalciumData(f);
        pupilData = io.getPupilData(f);
        noiseData = io.getVisNoiseInfo(f);
        rfData = io.getRFData(f);

        squW = diff(noiseData.edges(1:2)) / size(noiseData.frames,3);
        squH = diff(noiseData.edges(3:4)) / size(noiseData.frames,2);

        validPix = size(noiseData.frames, 3);
        if noiseData.edges(1) * noiseData.edges(2) < 0
            % determine right edge of all pixel columns
            rightEdges = noiseData.edges(1) + ...
                (1:size(noiseData.frames,3)) .* squW;
            validPix = find(rightEdges <= -60);
            noiseData.frames = noiseData.frames(:,:,validPix);
            noiseData.edges(2) = rightEdges(validPix(end));
        end

        % interpolate pupil data to match stimulus times
        t_stim = noiseData.times;
        eyePos = interp1(pupilData.time, pupilData.pos(:,1), t_stim);
        [~,~,bin] = histcounts(eyePos, ...
            prctile(eyePos, linspace(0,100,nb+1)));
        % get median eye position for each bin
        eyePosPerBin = NaN(1, nb);
        for b = 1:nb
            eyePosPerBin(b) = median(eyePos(bin == b), "omitnan");
        end

        % find neurons with good enough RFs to do shifts, ignore the rest
        validRF = find(rfData.pValues < minPVal(paramSet) & ...
            rfData.explVars > minExplainedVarianceStim(paramSet) & ...
            rfData.lambdas <= lambdas(paramSet));
        if length(validRF) < 1
            continue
        end
        % find neurons deemed good enough according to old criteria; fit
        % Gaussian to whole RF to determine median RF centre
        validOld = find(rfData.pValues < minPOld & ...
            rfData.explVars > minEVOld & ...
            rfData.lambdas <= lambdaOld);
        validOld = setdiff(validOld, validRF);

        %% Prepare calcium traces
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
                    nanInd2 = reshape(repmat(nanInd1, 1, upsample)', [], 1);
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
            y = traces(1:ind, indUnits(iUnit));
            y = fillmissing(y, 'linear');
            % fit double exponential to start of trace
            f = fit((1:length(y))', y, ...
                @(a,b,c,d,e,x) a + b .* exp(-x ./ c) + d .* exp(-x ./ e), ...
                'Lower', [0 0 0 0 0], ...
                'Upper', [max(y) max(y) 500 max(y) 500], ...
                'StartPoint', [min(y) mean(y) 50 mean(y) 5]);
            % remove fit
            traces(:, indUnits(iUnit)) = traces(:, indUnits(iUnit)) - ...
                f(1 : size(traces,1)) + f.a;
        end

        %% Map shifted RFs
        % map RFs for each eye position
        frames = noiseData.frames(noiseData.stimOrder,:,:);
        stimFrameDur = median(diff(t_stim));
        RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
            ceil(RFlimits(2) / stimFrameDur);
        RFs = NaN(length(validRF), size(frames,2), size(frames,3), length(RFtimesInFrames), ...
            2, nb);
        for p = 1:nb
            rfs = ...
                whiteNoise.getReceptiveField(traces, t_ca, frames, ...
                t_stim, RFtimesInFrames, lambdasStim(paramSet), crossFolds, ...
                bin ~= p);
            RFs(:,:,:,:,:,p) = permute(rfs, [5 1 2 3 4]);
        end
        % positive pixel in OFF subfield means: driven by black
        RFs(:,:,:,:,2,:) = -RFs(:,:,:,:,2,:);

        %% Fit RF maps with 2D Gaussian
        % perform gaussian fits on all newly mapped RFs
        % columns: amplitude, xSTD, yCentre, ySTD, rotation
        rfGaussConst = NaN(length(rfData.explVars), 5);
        % xCenter for each eye position
        rfGaussShift = NaN(length(rfData.explVars), nb);
        for iUnit = 1:length(validRF)
            % find best subfield (ON or OFF or ON/OFF), find best time of RF
            rf = squeeze(mean(RFs(iUnit,:,:,:,:,:),6)); % mean across shifts
            [~,mxTime] = max(max(reshape(permute(abs(rf), [1 2 4 3]), ...
                [], size(RFs,4)), [], 1));
            rf = squeeze(rf(:,:,mxTime,:));
            signs = NaN(1,2);
            subs = NaN(1,3);
            for sub = 1:2
                r = rf(:,:,sub);
                [subs(sub),ind] = max(abs(r), [], "all");
                signs(sub) = sign(r(ind));
            end
            rf(:,:,3) = (rf(:,:,1).*signs(1) + rf(:,:,2).*signs(2)) ./ 2;
            subs(3) = max(rf(:,:,3), [], "all");
            [m, mxSub] = max(subs);
            if subs(3) > 0.7 * m
                mxSub = 3;
            end
            if mxSub < 3
                rf = squeeze(RFs(iUnit,:,:,mxTime,mxSub,:)) .* signs(mxSub);
            else
                rf = squeeze(RFs(iUnit,:,:,mxTime,1,:) .* signs(1) + ...
                    RFs(iUnit,:,:,mxTime,2,:) .* signs(2)) ./ 2;
            end

            % interpolate RF so that pixels are square with edge length of 
            % 1 visual degree
            % note: xx and yy are locations of rf_visDeg relative to pixel
            % locations of rf (i.e. 1 is at the centre of the 1st pixel of
            % rf)
            [rf_visDeg,xx,yy] = whiteNoise.interpolateRFtoVisualDegrees(...
                rf, noiseData.edges);

            % fit 2D Gausian and shift parameter to RFs across eye positions
            fitPars = whiteNoise.fit2dGaussRF_wShift(rf_visDeg, false);
            % transform RF position to absolute values (relative to screen)
            fitPars(2) = noiseData.edges(1) + xx(1)*squW + fitPars(2);
            fitPars(4) = noiseData.edges(3) + yy(1)*squH + fitPars(4);

            % collect fitted parameters
            cellID = validRF(iUnit);
            rfGaussConst(cellID,:) = fitPars([1 3:6]);
            xCentres = fitPars(2) + [0 fitPars(7:end)];
            rfGaussShift(cellID,:) = xCentres;
        end

        %% Fit un-shifted RFs with 2D Gaussian 
        medianRFpositions = NaN(length(rfData.explVars), 1);
        % perform gaussian fit to determine x-position for old valid RFs
        for cellID = validOld'
            rf = squeeze(rfData.maps(cellID,:,validPix,:,:));
            [~,mxTime] = max(max(reshape(permute(abs(rf), [1 2 4 3]), ...
                [], size(rf,3)), [], 1));
            rf = squeeze(rf(:,:,mxTime,:));
            signs = NaN(1,2);
            subs = NaN(1,3);
            for sub = 1:2
                r = rf(:,:,sub);
                [subs(sub),ind] = max(abs(r), [], "all");
                signs(sub) = sign(r(ind));
            end
            rf(:,:,3) = (rf(:,:,1).*signs(1) + rf(:,:,2).*signs(2)) ./ 2;
            subs(3) = max(rf(:,:,3), [], "all");
            [m, mxSub] = max(subs);
            if subs(3) > 0.7 * m
                mxSub = 3;
                rf = rf(:,:,mxSub);
            else
                rf = rf(:,:,mxSub) .* signs(mxSub);
            end

            % interpolate RF so that pixels are square with edge length of 
            % 1 visual degree
            % note: xx and yy are locations of rf_visDeg relative to pixel
            % locations of rf (i.e. 1 is at the centre of the 1st pixel of
            % rf)
            [rf_visDeg,xx,yy] = whiteNoise.interpolateRFtoVisualDegrees(...
                rf, noiseData.edges);

            % fit 2D Gausian and shift parameter to RFs across eye positions
            fitPars = whiteNoise.fit2dGaussRF(rf_visDeg, false);
            % transform RF position to absolute values (relative to screen)
            fitPars(2) = noiseData.edges(1) + xx(1)*squW + fitPars(2);
            fitPars(4) = noiseData.edges(3) + yy(1)*squH + fitPars(4);

            rfGaussConst(cellID,:) = fitPars([1 3:6]);
            medianRFpositions(cellID) = fitPars(2);
        end
        
        %% Linear regression
        % linear regression relating pupil position in pixels to RF centre
        % in visual degrees &
        % determine median horizontal RF positions for newly fitted RFs,
        % only use units where all of the shifted x-positions are within
        % stimulus
        shifts = rfGaussShift; % [units x eyePosBins]
        shifts(shifts < noiseData.edges(1) | shifts > noiseData.edges(2)) = NaN;
        allShiftsValid = find(all(~isnan(shifts), 2));
        medianRFpositions(allShiftsValid) = median(shifts(allShiftsValid,:), 2);

        % fit: rf_xPos_norm = a + b * eye_xPos
        % rf_xPos_norm is different for each unit and is relative to its median
        % RF x-position (= rf_xPos - rf_medianXPos);
        % Vector eye_xPos is the same for each unit
        posNorm = shifts(allShiftsValid,:) - medianRFpositions(allShiftsValid);
        mdl = fitlm(reshape(repmat(eyePosPerBin, length(allShiftsValid), 1), [], 1), ...
            reshape(posNorm, [], 1), "RobustOpts", "on");
        coeffs = mdl.Coefficients.Estimate;

        % if some shifted RF positions were outside noise stimulus, use
        % linear model to estimate median RF position based on shifted RF
        % positions inside stimulus: 
        % rf_xPos - rf_medianXPos = a + b * eye_xPos
        % => rf_medianXPos = rf_xPos - a + b * eye_xPos
        someShiftsValid = setdiff(validRF, allShiftsValid);
        for cellID = someShiftsValid'
            x = shifts(cellID,:);
            ind = ~isnan(x);
            medianRFpositions(cellID) = mean(x(ind)' - predict(mdl, eyePosPerBin(ind)')); 
        end

        %% Save results
        % save RF fits
        folderRes = fullfile(folder.results, name, date);
        if ~isfolder(folderRes)
            mkdir(folderRes);
        end
        rfMapsShifted = NaN([size(rfData.maps,1), size(RFs,2:6)]);
        rfMapsShifted(validRF,:,:,:,:,:) = RFs;
        writeNPY(rfMapsShifted, fullfile(folderRes, 'rfPerEyePos.maps.npy'))
        writeNPY(rfGaussConst, fullfile(folderRes, 'rfPerEyePos.gaussFitPars_fixed.npy'))
        writeNPY(rfGaussShift, fullfile(folderRes, 'rfPerEyePos.gaussFitPars_xShifts.npy'))

        % save median horizontal RF positions
        writeNPY(medianRFpositions, fullfile(folderRes, 'rfPerEyePos.medianHorizRFPosition.npy'))
        isMeasured = false(length(caData.ids),1);
        isMeasured(allShiftsValid) = true;
        writeNPY(isMeasured, fullfile(folderRes, 'rfPerEyePos.isRFMedianMeasured.npy'))
        isModelled = false(length(caData.ids),1);
        isModelled(someShiftsValid) = true;
        writeNPY(isModelled, fullfile(folderRes, 'rfPerEyePos.isRFMedianModelled.npy'))
        isAcrossAllPupilPos = false(length(caData.ids),1);
        isAcrossAllPupilPos(validOld) = true;
        writeNPY(isAcrossAllPupilPos, fullfile(folderRes, 'rfPerEyePos.isRFMedianAcrossAllPupilPos.npy'))

        % save median eye positions
        writeNPY(eyePosPerBin, fullfile(folderRes, 'eyeToRFPos.eyeHorizPositions.npy'))

        % save linear regression fit
        writeNPY(coeffs, fullfile(folderRes, 'eyeToRFPos.regressionCoeffs.npy'))
        
        %% Make plots
        % plot RF maps and fitted ellipse for each neuron with shift
        for cellID = union(validRF, validOld)' %1:length(validRF)
            % find best subfield (ON or OFF or ON/OFF), find best time of RF
            if ismember(cellID, validRF)
                iUnit = find(validRF == cellID);
                rf = squeeze(mean(RFs(iUnit,:,:,:,:,:),6));
            else
                rf = squeeze(rfData.maps(cellID,:,validPix,:,:));
            end
            [~,mxTime] = max(max(reshape(permute(abs(rf), [1 2 4 3]), ...
                [], size(RFs,4)), [], 1));
            rf = squeeze(rf(:,:,mxTime,:));
            signs = NaN(1,2);
            subs = NaN(1,3);
            for sub = 1:2
                r = rf(:,:,sub);
                [subs(sub),ind] = max(abs(r), [], "all");
                signs(sub) = sign(r(ind));
            end
            rf(:,:,3) = (rf(:,:,1).*signs(1) + rf(:,:,2).*signs(2)) ./ 2;
            subs(3) = max(rf(:,:,3), [], "all");
            [m, mxSub] = max(subs);
            if subs(3) > 0.7 * m
                mxSub = 3;
            end

            if ismember(cellID, validRF)
                figure('Position', [4 300 1460 420])
                tiledlayout(2, nb, "TileSpacing", "tight")
                rf = squeeze(RFs(iUnit,:,:,mxTime,:,:));
                mx = max(abs(rf(:)));
                for sub = 1:2
                    for pos = 1:nb
                        nexttile
                        imagesc(...
                            [noiseData.edges(1)+squW/2 noiseData.edges(2)-squW/2], ...
                            [noiseData.edges(3)-squH/2 noiseData.edges(4)+squH/2], ...
                            rf(:,:,sub,pos),[-mx mx])
                        axis image off
                        colormap(gca, cms(:,:,sub))
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
                        n = x_rot < noiseData.edges(1) | x_rot > noiseData.edges(2);
                        x_rot(n) = NaN;
                        y_rot(n) = NaN;
                        plot(x_rot, y_rot, 'k')
                        if sub == 1
                            if pos == 1
                                title('nasal')
                            elseif pos == nb
                                title('temporal')
                            end
                        elseif sub == 2 && pos == 1
                            axis on
                            set(gca, 'box', 'off')
                        end
                    end
                    colorbar
                end
            else
                figure('Position', [400 300 700 420])
                tiledlayout(2, 1, "TileSpacing", "tight")
                rf = squeeze(rfData.maps(cellID,:,validPix,mxTime,:));
                mx = max(abs(rf(:)));
                for sub = 1:2
                    nexttile
                    imagesc(...
                        [noiseData.edges(1)+squW/2 noiseData.edges(2)-squW/2], ...
                        [noiseData.edges(3)-squH/2 noiseData.edges(4)+squH/2], ...
                        rf(:,:,sub),[-mx mx])
                    axis image off
                    colormap(gca, cms(:,:,sub))
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
                    n = x_rot < noiseData.edges(1) | x_rot > noiseData.edges(2) | ...
                        y_rot < noiseData.edges(3) | y_rot > noiseData.edges(4);
                    x_rot(n) = NaN;
                    y_rot(n) = NaN;
                    plot(x_rot, y_rot, 'k')
                    if sub == 2
                        axis on
                        set(gca, 'box', 'off')
                    end
                    colorbar
                end
            end
            sgtitle(sprintf('Neuron %d (plane %d, ID %d): %s (w/o shift: EV: %.3f, p: %.3f, lam: %.3f)', ...
                cellID, caData.planes(cellID), caData.ids(cellID), RFtypes{mxSub}, rfData.explVars(cellID), ...
                rfData.pValues(cellID), rfData.lambdas(cellID)), 'FontWeight', 'bold')
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
        hold on
        if ~isempty(allShiftsValid)
            plot(eyePosPerBin, median(shifts(allShiftsValid,:), 1, 'omitnan'), 'k', 'LineWidth', 2)
        end
        set(gca, 'box', 'off')
        axis tight
        xlabel('Pupil position (pixels)')
        ylabel('Horiz. RF centre (vis deg)')
        title('Median: included only if all shifts inside stimulus')

        if ~isempty(allShiftsValid)
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
