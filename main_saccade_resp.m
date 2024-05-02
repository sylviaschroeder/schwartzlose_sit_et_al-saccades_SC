function main_saccade_resp(folder, doPlot)

if nargin <2
    doPlot = 1;
end
% compute saccade responses and stats
%%
subjects = dir(fullfile(folder.data, 'SS*'));

%%
for subj = 4:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for iD = 1:length(dates)
        date = dates(iD).name;
 fprintf('Working on %s %s', name, date);
        %% load data
        this_folder = fullfile(folder.data, name, date, '001');
        gsData = io.getGrayScreenInfo(this_folder);

        if isempty(gsData.interval) % continue only if gray screen was recorded
            continue;
        end
   
        caData = io.getCalciumData(this_folder);
        valid_ca_frames = caData.time >= gsData.interval(1) & caData.time <= gsData.interval(2);
        caData.time = caData.time(valid_ca_frames);
        caData.traces = caData.traces(valid_ca_frames, :);

        pupilData = io.getPupilData(this_folder);
        valid_eye_frames = pupilData.time >= gsData.interval(1) & pupilData.time <= gsData.interval(2);
        pupilData.time = pupilData.time(valid_eye_frames);
        pupilData.pos = pupilData.pos(valid_eye_frames, :);
        pupilData.pupilSize = pupilData.pupilSize(valid_eye_frames, :);


        %% preproc ca data (align time across planes and detrend);
        neuralData = preprocCa(caData);
        nN = numel(caData.ids);

        %% convert eye movements to degrees
        this_folder  = fullfile(folder.results, name, date);
        try
            eyeModel = io.getEyeModel(this_folder);
        catch
            continue;
        end

        pupilData.degPosX = pupilData.pos(:,1)*eyeModel.coeff(2) + eyeModel.coeff(1);
        pupilData.degPosY = (pupilData.pos(:,2)-median(pupilData.pos(:,2), 'omitnan'))*eyeModel.coeff(2);

        % extract saccades, consider only movement on X
        %         [saccade_onoff, amplitudes, vel_stat, onsetXY] = ...
        %             eye.findSaccades(pupilData.pos(:,1), pupilData.pos(:,2), 3, 1, 'all',1);

        %% detect temp and nas saccade with same threshold (done in px
        % space, otherwise the fit changes for degrees cos bin size is
        % hardcoded in function
        % !!! add code step to clean up eye tracking from artefacts.

        [temp_saccade_onoff, temp_amplitudes, temp_vel_stat, temp_onsetXY] = ...
            eye.findSaccades(pupilData.degPosX, pupilData.degPosY, 3, 0.8, 'temp',1);

        [nas_saccade_onoff, nas_amplitudes, nas_vel_stat, nas_onsetXY] = ...
            eye.findSaccades(pupilData.degPosX, pupilData.degPosY, 3, 1, 'nas',1);

        %% combine nas and temp saccades 
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

        % get onsets and amplitude in degrees
        onset_deg = pupilData.degPosX(saccade_onoff(:,1));
        amplitude_deg = onset_deg + amplitudes.x'; % think about sign of coefficient

        %% build saccade matrix (currently only taking care of azimuth)
        %         RF_traj_X = pupilData.degPosX + eyeModel.medianRFpos';
        missingRF = isnan(eyeModel.medianRFpos(:,1));
        eyeModel.medianRFpos(missingRF,1) = median(eyeModel.medianRFpos(~missingRF,1), 'omitnan'); %remove when interpolated RF are provided
        onset_T = saccade_onoff(:,1);
        onset_RF = onset_deg + eyeModel.medianRFpos(:, 1)'; %nS*nN
        endpoint_RF = bsxfun(@plus, onset_RF,  amplitude_deg); %nS*nN

        saccade_matrix = zeros(numel(onset_deg), nN); % nSaccades * nN
        nas_valid = onset_RF> -125 & endpoint_RF < -55 & repmat(amplitude_deg, 1, nN)>0;
        temp_valid = onset_RF < -55 & endpoint_RF > -125  & repmat(amplitude_deg, 1, nN)<0;
        saccade_matrix(nas_valid) = 1;
        saccade_matrix(temp_valid) = -1;

        %% zscore using only moments outside of saccade response window
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

        %% saccade triggered averages
        [trial_nas, TA_nas, SE_nas, peak_nas, nas_peak_T, trial_window] = eye.saccade_ETA(neuralData, pupilData, onset_T, saccade_matrix, 'nas');
        [trial_temp, TA_temp, SE_temp, peak_temp, peak_temp_T] = eye.saccade_ETA(neuralData, pupilData, onset_T, saccade_matrix, 'temp');

        %% linear/circular shift test
        %compute linear shifts
        nShifts = 5000;
        [shifts, shift_onset_T, shift_lims] = eye.lin_shift(pupilData.time, onset_T, nShifts, 'circ');
        shift_saccade_matrix = saccade_matrix;
        shift_saccade_matrix(onset_T< shift_lims(1) | onset_T >= shift_lims(2), :) = [];

        % compute response statistic from permutations
        fprintf('Linear shifts...\n');
        % overfprintf(0, 'Working on shift n = ...')
        % nMsgChars =0;
        sh_peak_nas = zeros(nShifts, nN);
        sh_peak_temp = zeros(nShifts, nN);
        % fprintf('Progress:\n');
        % fprintf(['\n' repmat('.',1,nShifts) '\n\n']);
        parfor iSh = 1: nShifts
            [~, ~, ~, sh_peak_nas(iSh, :)] = eye.saccade_ETA(neuralData, pupilData, shift_onset_T(:, iSh), shift_saccade_matrix, 'nas');
            [~, ~, ~, sh_peak_temp(iSh, :)] = eye.saccade_ETA(neuralData, pupilData, shift_onset_T(:, iSh), shift_saccade_matrix, 'temp');
            % nMsgChars = overfprintf(nMsgChars, sprintf('%d', iSh))
            % fprintf('\b|\n');
        end
fprintf('Complete.\n');

        if subj == length(subjects) && iD==length(dates) % if last dataset, delete gcp
        delete(gcp('nocreate'));
        end

        % get a p-value from shifted statistic
        nas_p_up =[]; nas_p_dw =[];temp_p_up =[]; temp_p_dw =[];
        for iN = 1:nN
            nas_p_up(iN) = sum(sh_peak_nas(:, iN)>peak_nas(iN))/nShifts;
            nas_p_dw(iN) = sum(sh_peak_nas(:, iN)<peak_nas(iN))/nShifts;
            temp_p_up(iN) = sum(sh_peak_temp(:, iN)>peak_temp(iN))/nShifts;
            temp_p_dw(iN) = sum(sh_peak_temp(:, iN)<peak_temp(iN))/nShifts;

        end

        nas_p = min(nas_p_up, nas_p_dw)*2; % multiply by 2 because we run 2 test and taken the smallest
        missing = isnan(peak_nas);
        nas_p(missing) = NaN;

        temp_p = min(temp_p_up, temp_p_dw)*2;
        missing = isnan(peak_temp);
        temp_p(missing) = NaN;

        significant = nas_p<0.01 | temp_p <0.01;



        %% if requested, plot and save graphics
        if doPlot
         % deal with folders
        fprintf('  %s %s\n', name, date)
        folderPl = fullfile(folder.plots, 'Saccade_Responses', name, date);
        if ~isfolder(folderPl)
            mkdir(folderPl)
        end

         % plot the saccade matrix

            figure;
            eye.plot_saccade_matrix(saccade_matrix);
            print(fullfile(folder.plots, 'Saccade_Responses', name, date, sprintf('Saccade_matrix_%s_%s', ...
                name, date)), '-dpng'); close gcf
            % plot the raster of saccade responses across the population
            figure;
            eye.plot_saccade_raster(TA_nas, TA_temp, peak_nas, peak_temp);
            print(fullfile(folder.plots, 'Saccade_Responses', name, date, sprintf('All_neurons_STA_%s_%s', ...
                name, date)), '-dpng'); close gcf

            % plot the summary of the shift test
            figure;
            edges_p = 0:0.01:1;
            nas_p_density = histcounts(nas_p, edges_p);
            nas_p_density = nas_p_density/sum(nas_p_density);
            temp_p_density = histcounts(temp_p, edges_p);
            temp_p_density = temp_p_density/sum(temp_p_density);
            plot(edges_p(1:end-1), nas_p_density,'Linewidth', 2, 'Color', [0.8 0 0.3]); hold on
            plot(edges_p(1:end-1), temp_p_density,'Linewidth', 2, 'Color', [0 0.8 0.3]);
            xlim([-0.05 1.05]);
            ylim([0 0.15]);
            ylabel('Density')
            xlabel('shift test p-value')
            legend('Nas', 'Temp')
            formatAxes
            print(fullfile(folder.plots, 'Saccade_Responses', name, date, sprintf('Shift_test_values_%s_%s', ...
                name, date)), '-dpng'); close gcf

            % plot trial responses and shift test for significant neurons

            significant_idx = find(significant);
            figure;
            for iN = 1:numel(significant_idx)

                this_n = significant_idx(iN);
                subplot(2,2,1)
                if ~isempty(trial_nas{this_n})
                plot(trial_window, trial_nas{this_n }, 'Linewidth', 0.5, 'Color', [0.8 0 0.3 0.2]); hold on;
                plot(trial_window, TA_nas(:, this_n ), 'Linewidth', 2, 'Color', [0.8 0 0.3]);
                end
                ylabel('dF/F (zscore)')
                formatAxes

                subplot(2,2,2)
                sh_edges = [-2.5:0.1:2.5];
                sh_bins = sh_edges(1:end-1) +0.05;
                sh_density = histcounts(sh_peak_nas(:, this_n ), sh_edges);
                sh_density = sh_density/sum(sh_density);
                plot(sh_bins, sh_density, 'LineWidth',1, 'Color', [0.8 0 0.3]); hold on
                plot([peak_nas(this_n ), peak_nas(this_n )], [0 0.2], '--', 'LineWidth',1, 'Color', [0.8 0 0.3])
                ylim([0 0.25])
                ylabel('Density')
                title(sprintf('Nas p = %0.3f', nas_p(this_n)));
                formatAxes

                subplot(2,2,3)
                if ~isempty(trial_temp{this_n})
                plot(trial_window, trial_temp{this_n}, 'Linewidth', 0.5, 'Color', [0 0.8 0.3 0.2]); hold on;
                plot(trial_window, TA_temp(:, this_n), 'Linewidth', 2, 'Color', [0 0.8 0.3]);
                end
                ylabel('dF/F (zscore)')
                xlabel('time (s)')
                formatAxes

                subplot(2,2,4)
                sh_edges = [-2.5:0.1:2.5];
                sh_bins = sh_edges(1:end-1) +0.05;
                sh_density = histcounts(sh_peak_temp(:, this_n ), sh_edges);
                sh_density = sh_density/sum(sh_density);
                plot(sh_bins, sh_density, 'LineWidth',1, 'Color', [0 0.8 0.3]); hold on
                plot([peak_temp(this_n ), peak_temp(this_n )], [0 0.25], '--', 'LineWidth',1, 'Color', [0 0.8 0.3])
                ylim([0 0.25])
                ylabel('Density')
                xlabel('dF/F peak')
                title(sprintf('Temp p = %0.3f', temp_p(this_n)));
                formatAxes

                print(fullfile(folder.plots, 'Saccade_Responses', name, date, sprintf('Saccade_trials_neuron_%d', ...
                    this_n)), '-dpng');
                clf;
            end
            close gcf;

            % plot saccade TA for significant neurons
            figure;
            eye.plot_saccade_raster(TA_nas(:, significant), TA_temp(:, significant), peak_nas(significant), peak_temp(significant));
            print(fullfile(folder.plots, 'Saccade_Responses', name, date, sprintf('Significant_neurons_STA_%s_%s', ...
                name, date)), '-dpng'); close gcf
        end
        % % save results
        % folderRes = fullfile(folder.results, 'Saccade_Responses', name, date);
        % if ~isfolder(folderRes)
        %     mkdir(folderRes);
        % end
        %
        % writeNPY(onset_T, fullfile(folderRes, 'SaccadeResponses.onset_T.npy'))
        % writeNPY(onset_RF, fullfile(folderRes, 'SaccadeResponses.onset_RF.npy'))
        % writeNPY(endpoint_RF, fullfile(folderRes, 'SaccadeResponses.endpoint_RF.npy'))
        % writeNPY(saccade_matrix, fullfile(folderRes, 'SaccadeResponses.saccade_matrix.npy'))
        %
        % writeNPY(trial_nas, fullfile(folderRes, 'SaccadeResponses.trial_nas.npy'))
        % writeNPY(trial_temp, fullfile(folderRes, 'SaccadeResponses.trial_temp.npy'))
        % writeNPY(trial_window, fullfile(folderRes, 'SaccadeResponses.trial_window.npy'))
        %
        % writeNPY(peak_nas, fullfile(folderRes, 'SaccadeResponses.peak_nas.npy'))
        % writeNPY(peak_temp, fullfile(folderRes, 'SaccadeResponses.peak_temp.npy'))
        % writeNPY(sh_peak_nas, fullfile(folderRes, 'SaccadeResponses.sh_peak_nas.npy'))
        % writeNPY(sh_peak_temp, fullfile(folderRes, 'SaccadeResponses.sh_peak_temp.npy'))
        % writeNPY(nas_p, fullfile(folderRes, 'SaccadeResponses.nas_p.npy'))
        % writeNPY(temp_p, fullfile(folderRes, 'SaccadeResponses.temp_p.npy'))

    end
end

end