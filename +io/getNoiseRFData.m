function data = getNoiseRFData(folder)

data.maps = readNPY(fullfile(folder, '_ss_rf.maps.npy'));
data.explVars = readNPY(fullfile(folder, '_ss_rf.explVarsStim.npy'));
data.lambdas = readNPY(fullfile(folder, '_ss_rf.lambdasStim.npy'));
data.pValues = readNPY(fullfile(folder, '_ss_rf.pValues.npy'));
data.gaussPars = readNPY(fullfile(folder, '_ss_rf.gaussParameters.npy'));
data.peakToNoise = readNPY(fullfile(folder, '_ss_rf.peakToNoiseRatio.npy'));
data.bestSubfields = readNPY(fullfile(folder, '_ss_rf.bestSubField.npy'));
data.subfieldSigns = readNPY(fullfile(folder, '_ss_rf.subfieldSigns.npy'));
data.optimalDelays = readNPY(fullfile(folder, '_ss_rf.bestDelay.npy'));
data.timestamps = readNPY(fullfile(folder, '_ss_rfDescr.timestamps.npy'));
data.edges = readNPY(fullfile(folder, '_ss_rfDescr.edges.npy'));