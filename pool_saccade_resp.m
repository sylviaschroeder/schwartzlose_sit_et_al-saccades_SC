function pool_saccade_resp(folder, doPlot)

if nargin <2
    doPlot = 1;
end
fPlots = fullfile(folder.plots, 'saccadeVelocities');
if ~isfolder(fPlots)
    mkdir(fPlots)
end
% compute saccade responses and stats
%%
subjects = dir(fullfile(folder.data, 'SS*'));
% subjects = dir(fullfile(folder.data));
subjects = subjects(~startsWith({subjects.name},'.')); % remove . and  .. folders

%%
iSession = 0; 
for subj = 1:length(subjects)
% for subj = [1 9 10]
    name = subjects(subj).name;
    dates = dir(fullfile(folder.data, name, '2*'));
    for iD = 1:length(dates)
        iSession = iSession +1; 
        date = dates(iD).name;
        fprintf('Working on %s %s...\n', name, date);
        %% load data
        this_folder = fullfile(folder.results, name, date);
        [saccadeResponses(iSession), saccadeETA(iSession)] = io.getSaccadeResponses(folder);
         fprintf('...done \n', name, date);

    end
end

TA_nas = cat(1, saccadeETA(:).nasal);
TA_temp = cat(1, saccadeETA(:).temporal);

peak_nas = cat(1, saccadeResponses(:).peak_nas);
peak_temp = cat(1, saccadeResponses(:).peak_temporal);

  % figure;
  % eye.plot_saccade_raster(TA_nas(:, significant), TA_temp(:, significant), peak_nas(significant), peak_temp(significant), trial_window);
  % print(fullfile(folder.plots, 'Saccade_Responses', name, date, sprintf('Significant_neurons_STA_%s_%s', ...
  %     name, date)), '-dpng'); close gcf
    

end