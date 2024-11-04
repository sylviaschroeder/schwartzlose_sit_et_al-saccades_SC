function data = getRFFits(folder)

data.rfPosX = readNPY(fullfile(folder, 'rfPerEyePos.medianHorizRFPosition.npy'));
data.rfGaussPars = readNPY(fullfile(folder, 'rfPerEyePos.gaussFitPars_fixed.npy'));
data.isMeasured = readNPY(fullfile(folder, 'rfPerEyePos.isRFMedian_medianEyePos.npy'));
data.isModelled = readNPY(fullfile(folder, 'rfPerEyePos.isRFMedian_modelled.npy'));
if isfile(fullfile(folder, 'rfPerEyePos.gaussFitPars_xShifts.npy'))
    data.rfXShifts = readNPY(fullfile(folder, 'rfPerEyePos.gaussFitPars_xShifts.npy'));
else
    data.rfXShifts = [];
end