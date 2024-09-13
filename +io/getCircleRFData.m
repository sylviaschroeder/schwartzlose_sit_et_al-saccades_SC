function data = getCircleRFData(folder)

data.maps = readNPY(fullfile(folder, '_ss_circlesRf.maps.npy'));
data.explVars = readNPY(fullfile(folder, '_ss_circlesRf.explVarsStim.npy'));
data.lambdas = readNPY(fullfile(folder, '_ss_circlesRf.lambdasStim.npy'));
data.pValues = readNPY(fullfile(folder, '_ss_circlesRf.pValues.npy'));
data.timestamps = readNPY(fullfile(folder, '_ss_circlesRfDescr.timestamps.npy'));
data.x = readNPY(fullfile(folder, '_ss_circlesRfDescr.x.npy'));
data.y = readNPY(fullfile(folder, '_ss_circlesRfDescr.y.npy'));
data.diameters = readNPY(fullfile(folder, '_ss_circlesRfDescr.diameters.npy'));