function getRFposFromRetinotopy(folder)

%% Loop through datasets
folderPl = fullfile(folder.plots, 'RetinotopicRFs');
if ~isfolder(folderPl)
    mkdir(folderPl)
end

subjects = dir(fullfile(folder.data, 'SS*'));
for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        f = fullfile(folder.data, name, date, '001');
        fRes = fullfile(folder.results, name, date);

        if ~isfile(fullfile(f, '_ss_sparseNoise.times.npy'))
            continue
        end

        % load data
        brainPos = readNPY(fullfile(f, '_ss_2pRois.xyz.npy'));
        brainPos(:,3) = [];
        edges = readNPY(fullfile(f, '_ss_rfDescr.edges.npy'));
        edges(2) = min([edges(2) 0]);
        data = io.getRFFits(fRes);
        rf_x = data.rfPosX;
        rf_y = data.rfGaussPars(:,3);
        isMeasured = data.isMeasured;
        isModelled = data.isModelled;
        isAll = data.isAcrossAllPupilPos;
        
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
        folderRes = fullfile(folder.results, name, date);
        if ~isfolder(folderRes)
            mkdir(folderRes);
        end
        writeNPY([predict_rf_x, predict_rf_y], ...
            fullfile(folderRes, 'rfRetinotopy.pos.npy'))

        %% Make plots
        % figure('Position', [150 120 1640 770])
        figure('WindowState', 'maximized')
        tiledlayout(1,2)

        nexttile
        h = zeros(1,4);
        plot([rf_x(valid | isAll)'; predict_rf_x(valid | isAll)'], ...
            [rf_y(valid | isAll)'; predict_rf_y(valid | isAll)'], 'k')
        hold on
        h(1) = plot(rf_x(isAll), ...
            rf_y(isAll), ...
            'pentagram', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', [1 0.5 0.5], 'MarkerSize', 10);
        h(2) = plot(predict_rf_x(isAll), ...
            predict_rf_y(isAll), ...
            'o', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', [0.5 0.5 0.5]);
        h(3) = plot(rf_x(isMeasured | isModelled), ...
            rf_y(isMeasured | isModelled), ...
            'pentagram', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 10);
        h(4) = plot(predict_rf_x(isMeasured | isModelled), ...
            predict_rf_y(isMeasured | isModelled), ...
            'o', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', 'k');
        axis image
        set(gca, 'YDir', 'reverse', 'box', 'off')
        axis(edges)
        legend(h, 'Measured: bad RF', 'Fitted: bad RF', ...
            'Measured: good RF', 'Fitted: good RF')
        title('Measured vs fitted RF positions (only good RFs used for fit)')

        nexttile
        hold on
        h = [0 0];
        ind = isMeasured | isModelled | isAll;
        h(1) = plot(predict_rf_x(~ind), predict_rf_y(~ind), ...
            'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1].*0.75);
        h(2) = plot(predict_rf_x(ind), predict_rf_y(ind), ...
            'o', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 1 1].*0.25);
        axis image
        set(gca, 'YDir', 'reverse', 'box', 'off')
        legend(h, 'Unit without RF', 'Unit with RF')
        title('Fitted RF positions of all units')
        limits = [min(predict_rf_x) max(predict_rf_x) ...
            min(predict_rf_y) max(predict_rf_y)];
        limits([1 3]) = min([limits([1 3]); edges([1 3])], [], 1);
        limits([2 4]) = max([limits([2 4]); edges([2 4])], [], 1);
        axis(limits)
        if ~all(limits == edges)
            h(3) = plot(edges([1 1 2 2 1]), edges([3 4 4 3 3]), 'k:');
            legend(h, 'Unit without RF', 'Unit with RF', 'Monitor edge')
        end

        sgtitle(sprintf('%s %s', name, date), 'FontWeight', 'bold')
        saveas(gcf, fullfile(folderPl, sprintf('%s_%s.jpg', ...
            name, date)))
        close gcf
    end
end
