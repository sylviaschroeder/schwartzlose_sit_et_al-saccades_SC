function data = getPupilData(folder)

data.pupilSize = readNPY(fullfile(folder, 'eye.diameter.npy'));
data.time = readNPY(fullfile(folder, 'eye.timestamps.npy'));
data.pos = [];
if isfile(fullfile(folder, 'eye.xyPos.npy'))
    data.pos = readNPY(fullfile(folder, 'eye.xyPos.npy'));
end