function  [trial, TA, SE, peak, peak_T, window] = saccade_ETA(neuralData, pupilData, onset_T, saccade_matrix, type)

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
        this_neuron = neuralData.z_traces(:, iN);
        this_neuron(isnan(this_neuron)) = median(this_neuron, 'omitnan');
        this_saccades = saccade_matrix(:, iN);
        t_on = pupilData.time(onset_T(this_saccades))';
        [trial{iN}, TA(:, iN), SE(:, iN), window] = magicETA(neuralData.time, this_neuron, t_on, [-1 2], [-1 -0.5]);
        
        resp_win = window >resp_edges(1) & window <resp_edges(2);
        first = find(resp_win>0, 1, 'first');
        [~, peak_T(iN)] = max(abs(TA(resp_win, iN)));
        peak_T(iN) = peak_T(iN) + first - 1;
        peak(iN) = TA(peak_T(iN), iN);

    end

end