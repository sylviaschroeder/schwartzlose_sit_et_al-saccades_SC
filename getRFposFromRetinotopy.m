function getRFposFromRetinotopy(folder)

%% Loop through datasets
folderPl = fullfile(folder.plots, 'RetinotopicRFs');
if ~isfolder(folderPl)
    mkdir(folderPl)
end

subjects = dir(fullfile(folder.data));
subjects = subjects(~startsWith({subjects.name},'.'));

for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        f = fullfile(folder.data, name, date, '001');
        fRes = fullfile(folder.results, name, date);

        if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy')) && ...
                ~isfile(fullfile(f, "_ss_circles.intervals.npy"))
            continue
        end

        % load data
        brainPos = readNPY(fullfile(f, '_ss_2pRois.xyz.npy'));
        brainPos(:,3) = [];
        if isfile(fullfile(fRes, "_ss_circlesRf.pValues.npy"))
            rfData = io.getCircleRFData(fRes);
            gridW = median(diff(rfData.x));
            gridH = median(-diff(rfData.y));
            % edges: [left right top bottom] (above horizon: >0)
            edges = [rfData.x(1)-0.5*gridW rfData.x(end)+0.5*gridW ...
                rfData.y(1)+0.5*gridH rfData.y(end)-0.5*gridH];
            paradigm = 'circles';
        elseif isfile(fullfile(fRes, "_ss_rf.pValues.npy"))
            rfData = io.getNoiseRFData(fRes);
            % edges: [left right top bottom] (above horizon: >0)
            edges = rfData.edges;
            paradigm = 'visual noise';
        else
            continue
        end

        data = io.getRFFits(fRes);
        rf_x = data.rfPosX;
        rf_y = data.rfGaussPars(:,3);
        isMeasured = data.isMeasured;
        isModelled = data.isModelled;
        
        % for all units with a significant RF, use linear regression to
        % model:
        % rf_x = a1 + b1 * brain_x + c1 * brain_y
        % rf_y = a2 + b2 * brain_x + c2 * brain_y
        % valid = ~isnan(rf_x);
        valid = isMeasured | isModelled;
        fit_rf_x = fit(brainPos(valid,:), rf_x(valid), ...
            'poly11', 'Robust', 'Bisquare');
        fit_rf_y = fit(brainPos(valid,:), rf_y(valid), ...
            'poly11', 'Robust', 'Bisquare');
        % to directly check fitting, do those plots:
        % figure
        % plot(fit_rf_x, brainPos(valid,:), rf_x(valid))
        % xlabel('Brain in X')
        % ylabel('Brain in Y')
        % zlabel('RF azimuth')
        % title('Fit horiztonal RF position')
        % figure
        % plot(fit_rf_y, brainPos(valid,:), rf_y(valid))
        % xlabel('Brain in X')
        % ylabel('Brain in Y')
        % zlabel('RF elevation')
        % title('Fit vertical RF position')

        % for all units, use linear model to predict RF position
        predict_rf_x = fit_rf_x(brainPos);
        predict_rf_y = fit_rf_y(brainPos);

        %% Save results
        writeNPY([predict_rf_x, predict_rf_y], ...
            fullfile(fRes, 'rfRetinotopy.pos.npy'))

        %% Make plots
        % figure('Position', [150 120 1640 770])
        figure('WindowState', 'maximized')
        tiledlayout(1,2)

        nexttile
        h = zeros(1,2);
        plot([rf_x(valid)'; predict_rf_x(valid)'], ...
            [rf_y(valid)'; predict_rf_y(valid)'], 'k')
        hold on
        h(1) = plot(rf_x(valid), rf_y(valid), ...
            'pentagram', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 10);
        h(2) = plot(predict_rf_x(valid), predict_rf_y(valid), ...
            'o', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', 'k');
        axis image
        set(gca, 'YDir', 'normal', 'box', 'off')
        axis(edges([1 2 4 3]) + [-8.5 0 0 8.5])
        legend(h, 'Measured RF', 'Retinotopic RF')
        title('Measured vs retinotopic RF positions')

        nexttile
        hold on
        h = [0 0];
        h(1) = plot(predict_rf_x(~valid), predict_rf_y(~valid), ...
            'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1].*0.75);
        h(2) = plot(predict_rf_x(valid), predict_rf_y(valid), ...
            'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1].*0.25);
        axis image
        set(gca, 'YDir', 'normal', 'box', 'off')
        legend(h, 'Unit without RF', 'Unit with RF')
        title('Fitted RF positions of all units')
        limits = [min(predict_rf_x) max(predict_rf_x) ...
            min(predict_rf_y) max(predict_rf_y)];
        limits([1 3]) = min([limits([1 3]); edges([1 4])], [], 1);
        limits([2 4]) = max([limits([2 4]); edges([2 3])], [], 1);
        axis(limits)
        if ~all(limits == edges([1 2 4 3]))
            h(3) = plot(edges([1 1 2 2 1]), edges([3 4 4 3 3]), 'k:');
            legend(h, 'Unit without RF', 'Unit with RF', 'Monitor edge')
        end

        sgtitle(sprintf('%s %s (%s)', name, date, paradigm), ...
            'FontWeight', 'bold')
        saveas(gcf, fullfile(folderPl, sprintf('%s_%s.jpg', ...
            name, date)))
        close gcf
    end
end
