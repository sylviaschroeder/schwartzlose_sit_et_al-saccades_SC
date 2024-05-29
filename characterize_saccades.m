%visualize physical characteristics of saccades
function characterize_saccades(folder, doPlot)
if nargin <2
    doPlot = 1;
end
%%
subjects = dir(fullfile(folder.data, 'SS*'));
%%
for subj = 1:length(subjects)
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for iD = 1:length(dates)
        date = dates(iD).name;
        fprintf('Working on %s %s', name, date);
        folderPl = fullfile(folder.plots, 'Saccade_Characteristics', name, date);
        if ~isfolder(folderPl)
            mkdir(folderPl)
        end
        % load data
        this_folder = fullfile(folder.data, name, date, '001');

        gsData = io.getGrayScreenInfo(this_folder);

        if isempty(gsData.interval) % continue only if gray screen was recorded
            continue;
        end

        pupilData = io.getPupilData(this_folder);
        valid_eye_frames = pupilData.time >= gsData.interval(1) & pupilData.time <= gsData.interval(2);
        pupilData.time = pupilData.time(valid_eye_frames);
        pupilData.pos = pupilData.pos(valid_eye_frames, :);
        pupilData.pupilSize = pupilData.pupilSize(valid_eye_frames, :);

        % convert eye movements to degrees
        this_folder  = fullfile(folder.results, name, date);
        try
            eyeModel = io.getEyeModel(this_folder);
        catch
            continue;
        end

        pupilData.degPosX = pupilData.pos(:,1)*eyeModel.coeff(2) + eyeModel.coeff(1);
        pupilData.degPosY = (pupilData.pos(:,2)-median(pupilData.pos(:,2), 'omitnan'))*eyeModel.coeff(2);


        % extract all saccades
         [saccade_onoff, amplitudes, vel_stat, onsetXY] = ...
             eye.findSaccades(pupilData.degPosX, pupilData.degPosY, 3, 1, 'all',0);

        % detect temp and nas saccades
        [temp_saccade_onoff, temp_amplitudes, temp_vel_stat, temp_onsetXY] = ...
            eye.findSaccades(pupilData.degPosX, pupilData.degPosY, 3, 0.8, 'temp',0);

        [nas_saccade_onoff, nas_amplitudes, nas_vel_stat, nas_onsetXY] = ...
            eye.findSaccades(pupilData.degPosX, pupilData.degPosY, 3, 1, 'nas',0);

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

        % get onsets and amplitude in degrees
        onset_deg_temp_x = pupilData.degPosX(temp_saccade_onoff(:,1));
        onset_deg_temp_y = pupilData.degPosY(temp_saccade_onoff(:,2));
        amplitude_deg_temp_x = onset_deg_temp_x + temp_amplitudes.x';
        amplitude_deg_temp_y = onset_deg_temp_y + temp_amplitudes.y';

        onset_deg_nas_x = pupilData.degPosX(nas_saccade_onoff(:,1));
        onset_deg_nas_y = pupilData.degPosY(nas_saccade_onoff(:,2));
        amplitude_deg_nas_x = onset_deg_nas_x + nas_amplitudes.x';
        amplitude_deg_nas_y = onset_deg_nas_y + nas_amplitudes.y';

        onset_deg_x = pupilData.degPosX(saccade_onoff(:,1));
        onset_deg_y = pupilData.degPosY(saccade_onoff(:,2));
        amplitude_deg_x = onset_deg_x + amplitudes.x';
        amplitude_deg_y = onset_deg_y + amplitudes.y';

        %%characterize saccade physical dynamics
        % x = pupilData.pos(:,1);
        % y = pupilData.pos(:,2);

        %binning Y on X: amplitude depending on start position
        edges = min(onset_deg_x):5:max(onset_deg_x); %automate edge sizes
        [aveY, aveX, Xbin, stdY] = binYonX(onset_deg_x, amplitude_deg_x, edges);
        figure
        shadePlot(aveX, aveY, stdY, 'b')
        formatAxes()
        xlabel('start position in x')
        ylabel('x amplitude')
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Amplitude_vs_start_%s_%s', name, date)), '-dpng'); close gcf

        %same for nasal vs temporal
        startnasal = pupilData.degPosX(nas_saccade_onoff(:,1));
        startnasal_y = pupilData.degPosY(nas_saccade_onoff(:,1));
        starttemp = pupilData.degPosX(temp_saccade_onoff(:,1));
        starttemp_y = pupilData.degPosY(temp_saccade_onoff(:,1));
        endnasal = pupilData.degPosX(nas_saccade_onoff(:,2));
        endnasal_y = pupilData.degPosY(nas_saccade_onoff(:,2));
        endtemp = pupilData.degPosX(temp_saccade_onoff(:,2));
        endtemp_y = pupilData.degPosY(temp_saccade_onoff(:,2));
        edges = min(onset_deg_nas_x):5:max(onset_deg_nas_x);
        [aveY, aveX, Xbin, stdY] = binYonX(onset_deg_nas_x, amplitude_deg_nas_x, edges);
        figure
        shadePlot(aveX, aveY, stdY, 'green')
        hold on
        edges = min(onset_deg_temp_x):5:max(onset_deg_temp_x);
        [aveY, aveX, Xbin, stdY] = binYonX(onset_deg_temp_x, amplitude_deg_temp_x, edges);
        shadePlot(aveX, aveY, stdY, 'magenta')
        formatAxes()
        xlabel('start position in x')
        ylabel('x amplitude')
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Amplitude_vs_start_split_%s_%s', name, date)), '-dpng'); close gcf

        %findpeaks all
        [peakValues,peakLocations,widths,prominences] = findpeaks(vel_stat.velocity);
        figure
        scatter(widths,prominences)
        hold on
        xlabel('widths')
        ylabel('prominences')
        formatAxes()
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Find_peaks_all_%s_%s', name, date)), '-dpng'); close gcf

        %findpeaks nasal vs temporal
        [peakValues,peakLocations,widths,prominences] = findpeaks(nas_vel_stat.velocity);
        figure
        scatter(widths,prominences,'filled','green',MarkerFaceAlpha = 0.3)
        formatAxes()
        hold on
        [peakValues,peakLocations,widths,prominences] = findpeaks(temp_vel_stat.velocity);
        % %figure
        scatter(widths,prominences,'filled', 'magenta',MarkerFaceAlpha = 0.3)
        xlabel('widths')
        ylabel('prominences')
        formatAxes()
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Find_peaks_split_%s_%s', name, date)), '-dpng'); close gcf

        %all eye traces for nasal vs temporal
        etaT = pupilData.time(saccade_onoff(:,1));
        sample_rate = 1/30;
        min_time = -0.2;
        max_time = 0.4;
        periT = [min_time:sample_rate:max_time];
        [ETAmat_eye, ETA_eye, ETAse_eye, window_eye] = magicETA(pupilData.time, pupilData.degPosX, etaT, periT);
        etaT = pupilData.time(temp_saccade_onoff(:,1));
        [ETAmat_eye_temp, ETA_eye_temp, ETAse_eye_temp, window_eye_temp] = magicETA(pupilData.time, pupilData.degPosX, etaT, periT);
        ETA_eye_temp_mean = mean(ETA_eye_temp,2,'omitmissing');
        etaT = pupilData.time(nas_saccade_onoff(:,1));
        [ETAmat_eye_nas, ETA_eye_nas, ETAse_eye_nas, window_eye_nas] = magicETA(pupilData.time, pupilData.degPosX, etaT, periT);
        ETA_eye_nas_mean = mean(ETA_eye_nas,2,'omitmissing');
        figure
        subplot(2,1,1)
        plot(window_eye_temp,ETA_eye_temp,'Color',[1 0 1 0.3])
        hold on
        plot(window_eye_temp,ETA_eye_temp_mean,'m', 'LineWidth',3)
        title('temporal')
        ylabel('eye position in x')
        formatAxes()
        subplot(2,1,2)
        plot(window_eye_nas,ETA_eye_nas,'Color',[0 1 0 0.3])
        hold on
        plot(window_eye_nas,ETA_eye_nas_mean,'g', 'LineWidth',3)
        title('nasal')
        xlabel('time')
        ylabel('eye position in x')
        formatAxes()
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Traces_%s_%s', name, date)), '-dpng'); close gcf


        %histogram of amplitudes for nasal vs temporal
        figure
        histogram(abs(amplitude_deg_nas_x),'Facecolor','green')
        hold on
        histogram(abs(amplitude_deg_temp_x), 'Facecolor', 'magenta')
        xlabel('vector amplitude')
        ylabel('frequency')
        formatAxes()
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Histogram_Amplitudes_%s_%s', name, date)), '-dpng'); close gcf


        %define quadrants and angles
        difftemporalx = endtemp - starttemp;
        difftemporaly = endtemp_y - starttemp_y;
        diffnasalx = endnasal - startnasal;
        diffnasaly = endnasal_y - startnasal_y;
        nasal1 = find(diffnasaly>0);
        nasal1angle = atan(diffnasaly(nasal1)./(diffnasalx(nasal1)));
        nasal4 = find(diffnasaly<0);
        nasal4angle = atan(diffnasaly(nasal4)./(diffnasalx(nasal4)));
        temporal2 = find(difftemporaly>0);
        temporal2angle = pi + atan(difftemporaly(temporal2)./(difftemporalx(temporal2)));
        temporal3 = find(difftemporaly<0);
        temporal3angle = pi + atan(difftemporaly(temporal3)./(difftemporalx(temporal3)));
        nasalangles = cat(1,nasal1angle,nasal4angle);
        temporalangles = cat(1,temporal2angle,temporal3angle);

        %define peak speed for each saccade
        interval_intervals = [];
        peak_speed = [];
        for i = 1:height(saccade_onoff)
            interval_intervals = saccade_onoff(i,1):saccade_onoff(i,2);
            peak_speed(i) = max(vel_stat.velocity(interval_intervals(1,:)));
        end

        interval_intervals = [];
        nas_peak_speed = [];
        for i = 1:height(nas_saccade_onoff)
            interval_intervals = nas_saccade_onoff(i,1):nas_saccade_onoff(i,2);
            nas_peak_speed(i) = max(nas_vel_stat.velocity(interval_intervals(1,:)));
        end

        interval_intervals = [];
        temp_peak_speed = [];
        for i = 1:height(temp_saccade_onoff)
            interval_intervals = temp_saccade_onoff(i,1):temp_saccade_onoff(i,2);
            temp_peak_speed(i) = max(temp_vel_stat.velocity(interval_intervals(1,:)));
        end

        %binning Y on X: peak speed depending on start position
        edges = min(onset_deg_x):5:max(onset_deg_x);
        [aveY, aveX, Xbin, stdY] = binYonX(onset_deg_x, peak_speed', edges);
        figure
        shadePlot(aveX, aveY, stdY, 'b')
        formatAxes()
        xlabel('start position in x')
        ylabel('peak speed')
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('speed_vs_start__%s_%s', name, date)), '-dpng'); close gcf


        %histogram of speeds for nasal vs temporal
        edges = min(peak_speed):2:max(peak_speed);
        figure
        histogram(nas_peak_speed, edges, 'Facecolor', 'green')
        hold on
        histogram(temp_peak_speed, edges, 'Facecolor','magenta')
        xlabel('start position in x')
        ylabel('peak velocity')
        formatAxes()
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Histogram_velocity_%s_%s', name, date)), '-dpng'); close gcf


        %define the x and y components of the peak velocities for nasal vs temporal
        speedx_n = nas_peak_speed'.*(cos(nasalangles));
        speedy_n = nas_peak_speed'.*(sin(nasalangles));
        speedx_t = temp_peak_speed'.*(cos(temporalangles));
        speedy_t = temp_peak_speed'.*(sin(temporalangles));

        %plot quiver (subset to reduce overlap) and histogram of start position and velocity
        figure
        tiledlayout(2,1)
        ax1 = nexttile;
        quiver(onset_deg_nas_x(1:2:end),onset_deg_nas_y(1:2:end),speedx_n(1:2:end),speedy_n(1:2:end), 'off','g')
        hold on
        quiver(onset_deg_temp_x(1:2:end),onset_deg_temp_y(1:2:end),speedx_t(1:2:end),speedy_t(1:2:end), 'off','m')
        xlabel('position in x')
        ylabel('position in y')
        title('velocity')
        axis square
        formatAxes()
        ax2 = nexttile;
        edges = min(onset_deg_nas_x):5:max(onset_deg_nas_x);
        [aveY, aveX, Xbin, stdY] = binYonX(onset_deg_nas_x, nas_peak_speed', edges);
        shadePlot(aveX, aveY, stdY, 'g')
        hold on
        edges = min(onset_deg_temp_x):5:max(onset_deg_temp_x);
        [aveY, aveX, Xbin, stdY] = binYonX(onset_deg_temp_x, temp_peak_speed', edges);
        shadePlot(aveX, aveY, stdY, 'm')
        xlabel('starting position in x')
        ylabel('peak velocity')
        axis square
        formatAxes()
        linkaxes([ax1 ax2],'x')
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Quiver_velocity_%s_%s', name, date)), '-dpng'); close gcf


        %plot quiver (subset to reduce overlap) and histogram of start position and amplitude
        figure
        tiledlayout(2,1)
        ax1 = nexttile;
        quiver(onset_deg_nas_x(1:2:end),onset_deg_nas_y(1:2:end),amplitude_deg_nas_x(1:2:end),amplitude_deg_nas_y(1:2:end), 'off','g')
        hold on
        quiver(onset_deg_temp_x(1:2:end),onset_deg_temp_y(1:2:end),amplitude_deg_temp_x(1:2:end), amplitude_deg_temp_y(1:2:end), 'off','m')
        xlabel('position in x')
        ylabel('position in y')
        title('amplitude')
        axis square
        formatAxes()
        ax2 = nexttile;
        edges = min(onset_deg_nas_x):5:max(onset_deg_nas_x);
        [aveY, aveX, Xbin, stdY] = binYonX(onset_deg_nas_x, amplitude_deg_nas_x, edges);
        shadePlot(aveX, aveY, stdY, 'g')
        hold on
        edges = min(onset_deg_temp_x):5:max(onset_deg_temp_x);
        [aveY, aveX, Xbin, stdY] = binYonX(onset_deg_temp_x, amplitude_deg_temp_x, edges);
        shadePlot(aveX, aveY, stdY, 'm')
        xlabel('starting position in x')
        ylabel('amplitude in x')
        axis square
        formatAxes()
        linkaxes([ax1 ax2],'x')
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Quiver_amplitude_%s_%s', name, date)), '-dpng'); close gcf


        %compass plot of saccade endpoints relative to the same origin
        figure
        compass(amplitude_deg_temp_x,amplitude_deg_temp_y,'magenta')
        hold on
        compass(amplitude_deg_nas_x,amplitude_deg_nas_y,'green')
        formatAxes()
        title('x/y displacement')
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Compass_%s_%s', name, date)), '-dpng'); close gcf


        %polar histogram of saccade angles
        figure
        polarhistogram(nasalangles,'Facecolor','green')
        hold on
        polarhistogram(temporalangles,'Facecolor','magenta')
        formatAxes()
        title('frequency')
        print(fullfile(folder.plots, 'Saccade_Characteristics', name, date, sprintf('Polar_%s_%s', name, date)), '-dpng'); close gcf
        
        %save results
        folderRes = fullfile(folder.results, 'Saccade_Characteristics', name, date);
        if ~isfolder(folderRes)
            mkdir(folderRes);
        end

        % writeNPY(result_file, fullfile(folderRes, 'Saccade_Characteristics.results.npy'))

    end
end