
function data = getEyeModel(folder)


data.coeff = readNPY(fullfile(folder, 'eyeToRFPos.regressionCoeffs.npy'));
% data.medianRFpos = readNPY(fullfile(folder, 'rfPerEyePos.medianHorizRFPosition.npy'));
data.medianRFpos = readNPY(fullfile(folder, 'rfRetinotopy.pos.npy'));