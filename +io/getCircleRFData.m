function data = getCircleRFData(folder)

data.maps = readNPY(fullfile(folder, '_ss_circlesRf.maps.npy'));
data.explVars = readNPY(fullfile(folder, '_ss_circlesRf.explVarsStim.npy'));
data.lambdas = readNPY(fullfile(folder, '_ss_circlesRf.lambdasStim.npy'));
data.pValues = readNPY(fullfile(folder, '_ss_circlesRf.pValues.npy'));
data.gaussPars = readNPY(fullfile(folder, '_ss_circlesRf.gaussParameters.npy'));
data.peakToNoise = readNPY(fullfile(folder, '_ss_circlesRf.peakToNoiseRatio.npy'));
data.bestSubfields = readNPY(fullfile(folder, '_ss_circlesRf.bestSubField.npy'));
data.subfieldSigns = readNPY(fullfile(folder, '_ss_circlesRf.subfieldSigns.npy'));
data.optimalDelays = readNPY(fullfile(folder, '_ss_circlesRf.bestDelay.npy'));
data.sizeTuning = readNPY(fullfile(folder, '_ss_circlesRf.sizeTuning.npy'));
data.timestamps = readNPY(fullfile(folder, '_ss_circlesRfDescr.timestamps.npy'));
data.x = readNPY(fullfile(folder, '_ss_circlesRfDescr.x.npy'));
data.y = readNPY(fullfile(folder, '_ss_circlesRfDescr.y.npy'));
data.diameters = readNPY(fullfile(folder, '_ss_circlesRfDescr.diameters.npy'));