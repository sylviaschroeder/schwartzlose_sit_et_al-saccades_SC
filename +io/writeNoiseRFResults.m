function writeNoiseRFResults(results, folder)

writeNPY(results.maps, fullfile(folder, '_ss_rf.maps.npy'));
writeNPY(results.explVars, fullfile(folder, '_ss_rf.explVarsStim.npy'));
writeNPY(results.lambdas, fullfile(folder, '_ss_rf.lambdasStim.npy'));
writeNPY(results.pValues, fullfile(folder, '_ss_rf.pValues.npy'));
writeNPY(results.gaussPars, fullfile(folder, '_ss_rf.gaussParameters.npy'));
writeNPY(results.peakToNoise, fullfile(folder, '_ss_rf.peakToNoiseRatio.npy'));
writeNPY(results.bestSubfields, fullfile(folder, '_ss_rf.bestSubField.npy'));
writeNPY(results.subfieldSigns, fullfile(folder, '_ss_rf.subfieldSigns.npy'));
writeNPY(results.optimalDelays, fullfile(folder, '_ss_rf.bestDelay.npy'));
writeNPY(results.timestamps, fullfile(folder, '_ss_rfDescr.timestamps.npy'));
writeNPY(results.edges, fullfile(folder, '_ss_rfDescr.edges.npy'));