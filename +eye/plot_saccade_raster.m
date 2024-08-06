function plot_saccade_raster(nas_psth, temp_psth, nas_resp, temp_resp, trial_window)

% psth are nT by nN

clims = [-0.5 0.5];

nas_exc = abs(nas_resp) >= abs(temp_resp) & nas_resp >=0;
nas_inh = abs(nas_resp) >= abs(temp_resp) & nas_resp <0;

temp_exc = abs(nas_resp) < abs(temp_resp) & temp_resp >=0;
temp_inh = abs(nas_resp) < abs(temp_resp) & temp_resp <0;

[~, nas_sorting] = sort(nas_resp, 'descend');
[~, temp_sorting] = sort(nas_resp, 'descend');

ax(1) = subplot(4,2,1);
this_nas = nas_psth(:, nas_exc);
this_temp = temp_psth(:, nas_exc);
this_resp = nas_resp(nas_exc);
[sorted_nas, order] = sort_psth(this_nas, this_resp);
[sorted_temp, order] = sort_psth(this_temp, this_resp);

imagesc(trial_window, sum(nas_exc), sorted_nas'); colormap(colmaps.BlueWhiteRed);
caxis(clims)
ylabel('NA neurons')
title('Nas sacc');
colorbar;
formatAxes

ax(2) =subplot(4,2,2);
imagesc(trial_window, sum(nas_exc), sorted_temp'); colormap(colmaps.BlueWhiteRed);
caxis(clims)
title('Temp sacc');
colorbar;
formatAxes

ax(3) =subplot(4,2,3);

this_nas = nas_psth(:, nas_inh);
this_temp = temp_psth(:, nas_inh);
this_resp = nas_resp(nas_inh);
[sorted_nas, order] = sort_psth(this_nas, this_resp);
[sorted_temp, order] = sort_psth(this_temp, this_resp);

imagesc(trial_window, sum(nas_inh), sorted_nas'); colormap(colmaps.BlueWhiteRed);
caxis(clims)
formatAxes
ylabel('NS neurons')
colorbar;


ax(4) =subplot(4,2,4);
imagesc(trial_window, sum(nas_inh),sorted_temp'); colormap(colmaps.BlueWhiteRed);
caxis(clims)
colorbar;
formatAxes

ax(5) =subplot(4,2,5);

this_nas = nas_psth(:, temp_exc);
this_temp = temp_psth(:, temp_exc);
this_resp = temp_resp(temp_exc);
[sorted_nas, order] = sort_psth(this_nas, this_resp);
[sorted_temp, order] = sort_psth(this_temp, this_resp);

imagesc(trial_window, sum(temp_exc),sorted_nas'); colormap(colmaps.BlueWhiteRed);
caxis(clims)
formatAxes
ylabel('TA neurons')
colorbar;

ax(6) =subplot(4,2,6);
imagesc(trial_window, sum(temp_exc),sorted_temp'); colormap(colmaps.BlueWhiteRed);
caxis(clims)
colorbar;
formatAxes

ax(7) =subplot(4,2,7);

this_nas = nas_psth(:, temp_inh);
this_temp = temp_psth(:, temp_inh);
this_resp = temp_resp(temp_inh);
[sorted_nas, order] = sort_psth(this_nas, this_resp);
[sorted_temp, order] = sort_psth(this_temp, this_resp);

imagesc(trial_window, sum(temp_inh),sorted_nas'); colormap(colmaps.BlueWhiteRed);
caxis(clims)
formatAxes
xlabel('Time');
ylabel('TS neurons')
colorbar;


ax(8) =subplot(4,2,8);
imagesc(trial_window, sum(temp_inh),sorted_temp'); colormap(colmaps.BlueWhiteRed);
caxis(clims)
colorbar;
formatAxes
xlabel('Time');
colorbar;


linkaxes(ax(:), 'x');

end

function [sorted_psth, order] = sort_psth(psth, resp)

[~, order] = sort(resp, 'descend');

sorted_psth = psth(:, order);

end