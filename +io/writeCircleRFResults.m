function writeCircleRFResults(results, folder)

writeNPY(results.maps, fullfile(folder, '_ss_circlesRf.maps.npy'));
writeNPY(results.explVars, fullfile(folder, '_ss_circlesRf.explVarsStim.npy'));
writeNPY(results.lambdas, fullfile(folder, '_ss_circlesRf.lambdasStim.npy'));
writeNPY(results.pValues, fullfile(folder, '_ss_circlesRf.pValues.npy'));
writeNPY(results.gaussPars, fullfile(folder, '_ss_circlesRf.gaussParameters.npy'));
writeNPY(results.peakToNoise, fullfile(folder, '_ss_circlesRf.peakToNoiseRatio.npy'));
writeNPY(results.bestSubfields, fullfile(folder, '_ss_circlesRf.bestSubField.npy'));
writeNPY(results.subfieldSigns, fullfile(folder, '_ss_circlesRf.subfieldSigns.npy'));
writeNPY(results.optimalDelays, fullfile(folder, '_ss_circlesRf.bestDelay.npy'));
writeNPY(results.sizeTuning, fullfile(folder, '_ss_circlesRf.sizeTuning.npy'));
writeNPY(results.timestamps, fullfile(folder, '_ss_circlesRfDescr.timestamps.npy'));
writeNPY(results.x, fullfile(folder, '_ss_circlesRfDescr.x.npy'));
writeNPY(results.y, fullfile(folder, '_ss_circlesRfDescr.y.npy'));
writeNPY(results.diameters, fullfile(folder, '_ss_circlesRfDescr.diameters.npy'));