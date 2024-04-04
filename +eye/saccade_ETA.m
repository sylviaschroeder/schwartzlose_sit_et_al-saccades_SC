function  [trial, TA, SE, peak] = saccade_ETA(neuralData, pupilData, onset_T, saccade_matrix, type)

switch type
    case 'all'
        saccade_matrix = abs(saccade_matrix);
    case 'nas'
        saccade_matrix = saccade_matrix>0;
    case 'temp'
        saccade_matrix = saccade_matrix<0;

end


nN = numel(neuralData.ids);
    resp_edges = [-.5, 1.5];
    for iN = 1:nN
        t_on = pupilData.time(int16(onset_T(saccade_matrix(:, iN))))';
        [trial{iN}, TA(:, iN), SE(:, iN), window] = magicETA(neuralData.time, neuralData.traces(:,iN), t_on, [-1 3], [-1 -0.5]);
        
        resp_win = window >resp_edges(1) & window <resp_edges(2);
        peak(iN) = max(TA(resp_win, iN));
    end

end