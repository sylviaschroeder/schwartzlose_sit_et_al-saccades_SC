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
        fPlots = fullfile(folder.plots, 'RFs', name, date);
        if isfile(fullfile(f, "_ss_sparseNoise.times.npy"))
            fPlotsNoise = fullfile(folder.plots, 'RFs', name, date, ...
                'sorted_peakNoiseRatio_noise');
            if ~isfolder(fPlotsNoise)
                mkdir(fPlotsNoise);
            end
            results = io.getNoiseRFData(folderRes);
            valid = find(~isnan(results.peakToNoise));
            [~,sorted] = sort(results.peakToNoise(valid), "descend");
            sortedUnits = valid(sorted);
            for iUnit = 1:length(sortedUnits)
                file = sprintf('Unit%03d_noise.jpg', sortedUnits(iUnit));
                copyfile(fullfile(fPlots, file), ...
                    fullfile(fPlotsNoise, sprintf('%03d_%s', iUnit, file)));
            end
        end
        if isfile(fullfile(f, "_ss_circles.intervals.npy"))
            fPlotsCircles = fullfile(folder.plots, 'RFs', name, date, ...
                'sorted_peakNoiseRatio_circles');
            if ~isfolder(fPlotsCircles)
                mkdir(fPlotsCircles);
            end
            results = io.getCircleRFData(folderRes);
            valid = find(~isnan(results.peakToNoise));
            [~,sorted] = sort(results.peakToNoise(valid), "descend");
            sortedUnits = valid(sorted);
            for iUnit = 1:length(sortedUnits)
                file = sprintf('Unit%03d_circle.jpg', sortedUnits(iUnit));
                copyfile(fullfile(fPlots, file), ...
                    fullfile(fPlotsCircles, sprintf('%03d_%s', iUnit, file)));
            end
        end
    end
end
