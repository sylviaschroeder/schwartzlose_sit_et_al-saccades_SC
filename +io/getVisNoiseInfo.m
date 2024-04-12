function data = getGrayScreenInfo(folder)

data.stimOrder = readNPY(fullfile(folder, '_ss_sparseNoise._ss_sparseNoiseID.npy'));
data.times = readNPY(fullfile(folder, '_ss_sparseNoise.times.npy'));
data.edges = readNPY(fullfile(folder, '_ss_sparseNoiseArea.edges.npy'));
data.frames = readNPY(fullfile(folder, '_ss_sparseNoiseID.map.npy'));
data.interval = readNPY(fullfile(folder, '_ss_recordings.sparseNoise_intervals.npy'));