%% Folder definitions
folderData = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\DATA\DataToPublish\arousal_NYP_matlab\sc neurons 2p';
folderResults = 'D:\Sacccade Paper\Results';
folderCode = 'C:\dev';

%% Add paths
addpath(genpath(fullfile(folderCode, 'toolboxes\npy-matlab')))
addpath(genpath(fullfile(folderCode, 'workspaces\schwartzlose_sit_et_al-saccades_SC')))

%% Map eye position in pixels to eye position in visual degrees
main_eyePosInDegrees