function data = getCircleInfo(folder)

t = readNPY(fullfile(folder, '_ss_circles.intervals.npy'));
data.times = t(:,1);
data.diameter = readNPY(fullfile(folder, '_ss_circles.diameter.npy'));
data.isWhite = readNPY(fullfile(folder, '_ss_circles.isWhite.npy'));
data.xyPos = readNPY(fullfile(folder, '_ss_circles.xyPos.npy'));
data.interval = readNPY(fullfile(folder, '_ss_recordings.circles_intervals.npy'));