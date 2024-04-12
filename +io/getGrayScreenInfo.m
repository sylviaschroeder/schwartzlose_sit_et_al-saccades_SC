function data = getGrayScreenInfo(folder)

if exist(fullfile(folder, '_ss_recordings.grayScreen_intervals.npy'), 'file')

    data.interval = readNPY(fullfile(folder, '_ss_recordings.grayScreen_intervals.npy'));
else

    data.interval =[];
    warning('No gray screen was recorded');

end

