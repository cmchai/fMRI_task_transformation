% Script to get data(contrasts from individual subj) from ROIs for Task transformation paradigm
% Author: Mengqiao Chai (Mengqiao.Chai@ugent.be)
% Last update on: 9th, Jan, 2024

clear all; close all;
addpath(genpath('/Users/mengqiao/Documents/fMRI_task_transform/MRI_scripts')); % to add all scripts, incl toolboxes + spm12
addpath(genpath('/Users/mengqiao/Documents/MATLAB/packages/spm12'));

% ROI directory
ROI_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/ROI/AAL3';
ROI_mat = 'AAL3_rois.mat';
ROIs = who('-file', fullfile(ROI_dir, ROI_mat)); % return a cell array of ROI variable names
load(fullfile(ROI_dir, ROI_mat))

% data directories
SLM_root_dir = '/Users/mengqiao/Documents/fMRI_task_transform/MRI_data/Task_transform/second_level/univariate';
Model = 'FIR-01B';
Model_SLM_root_dir = fullfile(SLM_root_dir, Model, 'results');

% contrasts to look at in the data directory
Contrasts = [21:40];

% define number of participants
sum_part = 43;

% create data matrix
data_array = zeros(sum_part, length(Contrasts), length(ROIs)); % number of participants X number of contrasts X number of ROIs

for i = 1:length(ROIs) % loop over each ROI

    currentROI = evalin('base', ROIs{i});

    for ii = 1:length(Contrasts)

        contrast = ['Contrast', num2str(Contrasts(ii))];
        contrast_dir = fullfile(Model_SLM_root_dir, contrast);

        % load SPM mat
        clear SPM    % clear the SPM file from the previous participant   
        load(fullfile(contrast_dir, 'SPM.mat'));
        data_mean = mean(spm_get_data(SPM.xY.P, currentROI),2); % mean across voxels of a certain ROI
        data_array(:,ii,i) = data_mean;
    end
end


data_array_mean = squeeze(mean(data_array, 1));
data_array_sd = squeeze(std(data_array,0,1));

x = [1:10];
y_rg = data_array_mean(1:10,1);
sd_rg = data_array_sd(1:10,1);

y_tr = data_array_mean(11:20,1);
sd_tr = data_array_sd(11:20,1);


figure
plot(x, y_rg, 'b-o', x, y_tr, 'r-o', "LineWidth", 1.5)
yline(0, '--')


figure
errorbar(x, y_rg, sd_rg, 'b-o', "LineWidth", 1.5, "MarkerSize",12,...
    "MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
hold on
errorbar(x+0.05, y_tr, sd_tr, 'r-o', "LineWidth", 1.5, "MarkerSize",12,...
    "MarkerEdgeColor","red","MarkerFaceColor",[0.8500, 0.3250, 0.0980])
ylim([-0.03 0.03])
yline(0, '--')


