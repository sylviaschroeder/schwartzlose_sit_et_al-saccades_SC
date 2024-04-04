
% compute saccade responses and stats
%%
subjects = dir(fullfile(folder.data, 'SS*'));

for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for dt = 1:length(dates)
        date = dates(dt).name;
        this_folder = fullfile(folder.data, name, date, '001');
        
        % load data
        caData = io.getCalciumData(this_folder);
        pupilData = io.getPupilData(this_folder);

        % preproc ca data (align time across planes and detrend);
        neuralData = preprocCa(caData);
        nN = numel(caData.ids);

        % convert eye movements to degrees
        this_folder  = fullfile(folder.results, 'RFsByEyePos', name, date);
        eyeModel = io.getEyeModel(this_folder);
      
        pupilData.degPosX = pupilData.pos(:,1)*eyeModel.coeff(2) + eyeModel.coeff(1);
        pupilData.degPosY = zeros(size(pupilData.degPosX));

        % extract saccades, consider only movement on X
        [saccade_onoff, amplitudes, vel_stat, onsetXY] = ...
            eye.findSaccades(pupilData.pos(:,1), pupilData.pos(:,2), 3, 1, 'all',1);

%         [saccade_onoff, amplitudes, vel_stat, onsetXY] = ...
%             eye.findSaccades(pupilData.pos(:,1), pupilData.pos(:,2), 3, 1, 'temp',1);
% 
%         [saccade_onoff, amplitudes, vel_stat, onsetXY] = ...
%             eye.findSaccades(pupilData.pos(:,1), pupilData.pos(:,2), 3, 1, 'nas',1);

        % convert eye movements to degrees
        onset_deg = pupilData.degPosX(saccade_onoff(:,1));
        amplitude_deg = onset_deg + amplitudes.x'*eyeModel.coeff(2);

        % build saccade matrix
        %         RF_traj_X = pupilData.degPosX + eyeModel.medianRFpos';
        onset_T = saccade_onoff(:,1);
        onset_RF = onset_deg + eyeModel.medianRFpos';
        endpoint_RF = bsxfun(@plus, onset_RF,  amplitude_deg);

        saccade_matrix = zeros(numel(onset_deg), nN); % nSaccades * nN
        nas_valid = onset_RF> -125 & endpoint_RF < -55 & amplitude_deg>0;
        temp_valid = onset_RF < -55 & endpoint_RF > -125  & amplitude_deg<0;
        saccade_matrix(nas_valid) = 1; 
        saccade_matrix(temp_valid) = -1;
        figure; imagesc(saccade_matrix); 

        % saccade triggered averages
        [nas_trial, nas_TA, nas_SE, nas_peak] = eye.saccade_ETA(neuralData, pupilData, onset_T, saccade_matrix, 'nas');
        [temp_trial, temp_TA, temp_SE, temp_peak] = eye.saccade_ETA(neuralData, pupilData, onset_T, saccade_matrix, 'temp');

        % compute linear shifts
       
        nShifts = 10;
        nTs = numel(pupilData.time);
        shift_lims = round(prctile(onset_T(:,1), [10 90]));

        shift_neg = rand(round(nShifts/2), 1)*(-1)*shift_lims(1) +1;
        shift_pos = rand(round(nShifts/2), 1)*(nTs - shift_lims(2)) - 1;
        shifts = [0; shift_neg; shift_pos];
        
        shift_onset_T = onset_T;
        shift_onset_T(onset_T< shift_lims(1) | onset_T >= shift_lims(2)) =[];
        shift_onset_T = shift_onset_T + shifts'; % nSaccade * nShifts
        shift_saccade_matrix = saccade_matrix;
        shift_saccade_matrix(onset_T< shift_lims(1) | onset_T >= shift_lims(2), :) = [];

        % linear shift test
        
        for iSh = 1: nShifts
            [~, ~, ~, nas_peak(iSh, :)] = eye.saccade_ETA(neuralData, pupilData, shift_onset_T(:, iSh), shift_saccade_matrix, 'nas');
            [~, ~, ~, temp_peak(iSh, :)] = eye.saccade_ETA(neuralData, pupilData, shift_onset_T(:, iSh), shift_saccade_matrix, 'temp');

        end


    end
end

 