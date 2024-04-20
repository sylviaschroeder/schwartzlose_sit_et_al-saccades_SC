function data = getRFFits(folder)

data.rfPosX = readNPY(fullfile(folder, 'rfPerEyePos.medianHorizRFPosition.npy'));
data.rfGaussPars = readNPY(fullfile(folder, 'rfPerEyePos.gaussFitPars_fixed.npy'));
data.rfXShifts = readNPY(fullfile(folder, 'rfPerEyePos.gaussFitPars_xShifts.npy'));
data.isMeasured = readNPY(fullfile(folder, 'rfPerEyePos.isRFMedianMeasured.npy'));
data.isModelled = readNPY(fullfile(folder, 'rfPerEyePos.isRFMedianModelled.npy'));
data.isAcrossAllPupilPos = readNPY(fullfile(folder, 'rfPerEyePos.isRFMedianAcrossAllPupilPos.npy'));