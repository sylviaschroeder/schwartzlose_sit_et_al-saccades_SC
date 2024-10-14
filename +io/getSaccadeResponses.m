function [saccadeResponse, saccadeETA] = getSaccadeResponses(folder)

if exist(fullfile(folder, 'saccadeResponses.intervals.npy'), 'file')

    % about saccades: [nEvents x ...]
    % saccade.intervals
    saccadeResponse.intervals = readNPY(fullfile(folder, 'saccadeResponses.intervals.npy')); % nEvents x 2
    % saccade.amplitude
    saccadeResponse.amplitudes = readNPY(fullfile(folder, 'saccadeResponses.saccadeAmplitude.npy')); % nEvents x 2

    % about neurons: [nNeurons x ...]
    % saccadeResponse.RFPosAtStart (you mean the position of the RF at saccade start?)
    saccadeResponse.RFPosAtStart = readNPY(fullfile(folder, 'saccadeResponse.RFPosAtStart.npy')); %nEvents x nNeurons
    % saccadeResponse.RFPosAtEnd
    saccadeResponse.RFPosAtEnd = readNPY(fullfile(folder, 'saccadeResponse.RFPosAtEnd.npy')); %nEvents x nNeurons
    % saccadeResponse.valid
    saccadeResponse.valid = readNPY(fullfile(folder, 'saccadeResponse.valid.npy')); %nEvents x nNeurons
    % saccadeResponse.peakNasal
    saccadeResponse.peakNasal = readNPY(fullfile(folder, 'saccadeResponse.peakNasal.npy')); % 1x nNeurons
    % saccadeResponse.peakTemporal
    saccadeResponse.peakTemporal = readNPY(fullfile(folder, 'saccadeResponse.peakTemporal.npy')) ; % 1x nNeurons
    % saccadeResponse.peakNasalShuffled
    saccadeResponse.peakNasalShuffled = readNPY(fullfile(folder, 'saccadeResponse.peakNasalShuffled.npy'));  % nShuffles x nNeurons
    % saccadeResponse.peakTemporalShuffled
    saccadeResponse.peakTemporalShuffled = readNPY(fullfile(folder, 'saccadeResponse.peakTemporalShuffled.npy')); % nShuffles x nNeurons
    % saccadeResponse.pValueNasal
    saccadeResponse.pValueNasal = readNPY(fullfile(folder, 'saccadeResponse.pValueNasal.npy')); % 1 x nNeurons
    % saccadeResponse.pValueTemporal
    saccadeResponse.pValueTemporal = readNPY(fullfile(folder, 'saccadeResponse.pValueTemporal.npy')); % 1 x nNeurons

    % about saccade response ETA/kernel: [t x ...]
    % saccadeETA.nasal
    saccadeETA.nasal = readNPY(fullfile(folder, 'saccadeETA.nasal.npy')); % trialWindowT x nNeurons
    % saccadeETA.temporal
    saccadeETA.temporal = readNPY(fullfile(folder, 'saccadeETA.temporal.npy')); % trialWindowT x nNeurons
    % saccadeETA.timestamps
    saccadeETA.timestamps = readNPY(fullfile(folder, 'saccadeETA.timestamps.npy')); %nT x 1;

else

    saccadeResponse.intervals = [];
    saccadeResponse.amplitudes =  [];
    saccadeResponse.RFPosAtStart =  [];
    saccadeResponse.RFPosAtEnd =  [];
    saccadeResponse.valid =  [];
    saccadeResponse.peakNasal =  [];
    saccadeResponse.peakTemporal =  [];
    saccadeResponse.peakNasalShuffled =  [];
    saccadeResponse.peakTemporalShuffled =  [];
    saccadeResponse.pValueNasal =  [];
    saccadeResponse.pValueTemporal =  [];
    saccadeETA.nasal =  [];
    saccadeETA.temporal =  [];
    saccadeETA.timestamps =  [];

    warning('no saccade responses found');

end
