function data = getCalciumData(folder)

data.planes = readNPY(fullfile(folder, '_ss_2pRois._ss_2pPlanes.npy'));
data.ids = readNPY(fullfile(folder, '_ss_2pRois.ids.npy'));
data.traces = readNPY(fullfile(folder, '_ss_2pCalcium.dff.npy'));
data.time = readNPY(fullfile(folder, '_ss_2pCalcium.timestamps.npy'));
data.delays = readNPY(fullfile(folder, '_ss_2pPlanes.delay.npy'));