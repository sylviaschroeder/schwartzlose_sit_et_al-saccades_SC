%% Folder definitions
% folders Sylvia (please comment out and add your own lines)
folder.data = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\DATA\DataToPublish\arousal_NYP_matlab\sc neurons 2p';
folder.results = 'D:\Sacccade Paper\Results';
folder.plots = 'D:\Sacccade Paper\Plots';
folder.codeToolboxes = 'C:\dev\toolboxes';
folder.codeThis = 'C:\dev\workspaces';

%% Add paths
addpath(genpath(fullfile(folder.codeToolboxes, 'npy-matlab')))
addpath(genpath(fullfile(folder.codeThis, 'schwartzlose_sit_et_al-saccades_SC')))

%% Map eye position in pixels to eye position in visual degrees
determineEyePosInDegrees(folder);