
% compute saccade responses and stats
%%
subjects = dir(fullfile(folder.data, 'SS*'));

for subj = 2%1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for iD = 1:length(dates)
        date = dates(iD).name;
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
%         [saccade_onoff, amplitudes, vel_stat, onsetXY] = ...
%             eye.findSaccades(pupilData.pos(:,1), pupilData.pos(:,2), 3, 1, 'all',1);

        % detect temp and nas saccade with same threshold (done in px
        % space, otherwise the fit changes for degrees cos bin size is
        % hardcoded in function

        [temp_saccade_onoff, temp_amplitudes, temp_vel_stat, temp_onsetXY] = ...
            eye.findSaccades(pupilData.pos(:,1), pupilData.pos(:,2), 3, 1, 'temp',1);

        [nas_saccade_onoff, nas_amplitudes, nas_vel_stat, nas_onsetXY] = ...
            eye.findSaccades(pupilData.pos(:,1), pupilData.pos(:,2), 3, 1, 'nas',1);

        % combine saccades 
        saccade_onoff = cat(1, temp_saccade_onoff, nas_saccade_onoff);
        onsetXY = cat(1, temp_onsetXY, nas_onsetXY);
        amplitudes.vec = cat(2, temp_amplitudes.vec, nas_amplitudes.vec);
        amplitudes.x = cat(2, temp_amplitudes.x, nas_amplitudes.x);
        amplitudes.y = cat(2, temp_amplitudes.y, nas_amplitudes.y);

        [~, sortidx] = sort(saccade_onoff(:,1), 'ascend');
        saccade_onoff = saccade_onoff(sortidx,:);
        onsetXY = onsetXY(sortidx,:);
        amplitudes.vec = amplitudes.vec(sortidx);
        amplitudes.x = amplitudes.x(sortidx);
        amplitudes.y= amplitudes.y(sortidx);

        % convert eye movements to degrees
        onset_deg = pupilData.degPosX(saccade_onoff(:,1));
        amplitude_deg = onset_deg + amplitudes.x'*eyeModel.coeff(2); % think about sign of coefficient

        % build saccade matrix
        %         RF_traj_X = pupilData.degPosX + eyeModel.medianRFpos';
        missingRF = isnan(eyeModel.medianRFpos);
        eyeModel.medianRFpos(missingRF) = median(eyeModel.medianRFpos(~missingRF), 'omitnan'); %remove when interpolated RF are provided
        onset_T = saccade_onoff(:,1);
        onset_RF = onset_deg + eyeModel.medianRFpos';
        endpoint_RF = bsxfun(@plus, onset_RF,  amplitude_deg);

        saccade_matrix = zeros(numel(onset_deg), nN); % nSaccades * nN
        nas_valid = onset_RF> -1225 & endpoint_RF < -5 & amplitude_deg>0;
        temp_valid = onset_RF < -5 & endpoint_RF > -1225  & amplitude_deg<0;
        saccade_matrix(nas_valid) = 1; 
        saccade_matrix(temp_valid) = -1;
        figure; imagesc(saccade_matrix); 

        % zscore neural activity;

        dT = median(diff(neuralData.time)); % frame rate
        frame_resp_window = round(-0.5/dT):round(2/dT); % index of frames around event times
        saccade_on_in_frames =[];
       for iS = 1:numel(onset_T)
        [~, saccade_on_in_frames(iS)] = min(abs(neuralData.time - pupilData.time(onset_T(iS)))); 
       end
        saccade_on_in_frames = saccade_on_in_frames + frame_resp_window';
        saccade_off_in_frames =setdiff(1:numel(neuralData.time), saccade_on_in_frames(:));
        for iN = 1:nN
            neuralData.z_traces(:, iN) = (neuralData.traces(:, iN) - ...
                nanmean(neuralData.traces(saccade_off_in_frames(:), iN))) /...
                nanstd(neuralData.traces(saccade_off_in_frames(:), iN));
        end

        % saccade triggered averages
        [nas_trial, nas_TA, nas_SE, nas_peak, nas_peak_T] = eye.saccade_ETA(neuralData, pupilData, onset_T, saccade_matrix, 'nas');
        [temp_trial, temp_TA, temp_SE, temp_peak, temp_peak_T] = eye.saccade_ETA(neuralData, pupilData, onset_T, saccade_matrix, 'temp');

        [~, sort_nas] = sort(nas_peak, 'descend');
        [~, sort_temp] = sort(temp_peak, 'descend');

        figure; 
        subplot(1,2,1)
        imagesc(nas_TA(:, sort_nas)');
        caxis([-0.2 0.2])
        subplot(1,2,2)
        imagesc(temp_TA(:, sort_temp)');
        caxis([-0.2 0.2])

        % compute linear shifts
       
        nShifts = 1000;
        nTs = numel(pupilData.time);
        shift_lims = round(prctile(onset_T(:,1), [10 90]));

        shift_neg = floor(rand(round(nShifts/2), 1)*(-1)*shift_lims(1)) +1; % make shift proportional
        shift_pos = floor(rand(round(nShifts/2), 1)*(nTs - shift_lims(2))) - 1;
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

 