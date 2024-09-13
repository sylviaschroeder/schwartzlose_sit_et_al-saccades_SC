function writeNoiseRFResults(results, folder)

writeNPY(results.maps, fullfile(folder, '_ss_rf.maps.npy'));
writeNPY(results.explVars, fullfile(folder, '_ss_rf.explVarsStim.npy'));
writeNPY(results.lambdas, fullfile(folder, '_ss_rf.lambdasStim.npy'));
writeNPY(results.pValues, fullfile(folder, '_ss_rf.pValues.npy'));
writeNPY(results.timestamps, fullfile(folder, '_ss_rfDescr.timestamps.npy'));
writeNPY(results.edges, fullfile(folder, '_ss_rfDescr.edges.npy'));