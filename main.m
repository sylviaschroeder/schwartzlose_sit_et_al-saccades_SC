
clear; 

%% Folder definitions
% folders Sylvia (please comment out and add your own lines)
% folder.data = 'C:\Users\Sylvia\OneDrive - University of Sussex\Lab\DATA\DataToPublish\arousal_NYP_matlab\sc neurons 2p';
% folder.results = 'D:\Sacccade Paper\Results';
% folder.plots = 'D:\Sacccade Paper\Plots';
% folder.codeToolboxes = 'C:\dev\toolboxes';
% folder.codeThis = 'C:\dev\workspaces';

% % folders Federico Laptop
% folder.data = '/Users/federico/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Neural correlates of eye movements/Paper/SC_data';
% folder.results = '/Users/federico/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Neural correlates of eye movements/Paper/Results';
% folder.plots = '/Users/federico/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Neural correlates of eye movements/Paper/Plots';
% folder.codeToolboxes = '/Users/federico/Documents/GitHub';
% folder.codeThis = '/Users/federico/Documents/GitHub';
% 

% folders Federico's WS
folder.data = 'D:\OneDrive - University College London\Neural correlates of eye movements\Paper\SC_data';
folder.results = 'D:\OneDrive - University College London\Neural correlates of eye movements\Paper\Results';
folder.plots = 'D:\OneDrive - University College London\Neural correlates of eye movements\Paper\Plots';
folder.codeToolboxes = 'D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\Stable';
folder.codeThis = 'D:\OneDrive - Fondazione Istituto Italiano Tecnologia\Documents\Code\Dev';


%% Add paths
addpath(genpath(fullfile(folder.codeToolboxes, 'npy-matlab')))
addpath(genpath(fullfile(folder.codeThis, 'schwartzlose_sit_et_al-saccades_SC')))
addpath(genpath(fullfile(folder.codeToolboxes, 'FedBox')))

%% Map eye position in pixels to eye position in visual degrees (Sylvia)
determineEyePosInDegrees(folder);

%% Saccade dynamics (Heather)


%% Saccade triggered responses (Federico)
main_saccade_resp(folder);

%% REVERT TEST