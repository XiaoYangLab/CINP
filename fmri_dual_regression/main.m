%% generate the functional activity amplitude and functional connectivity of RSNs using dual regression
% Author: Yang Xiao, PKU, 2024
% @: xiaoyang9604@gmail.com

close all; 
clear all;
clc;

%% set working place
data_idx = 'ZMD'; % PKU6\\\ZMD
data_dir = 'F:\fAPCS\workpipeline\SZImagingGenetic\0_dataset\Output\APDs_fMRI\';
data_name = ['SZ_', data_idx, '_rsfMRI_baseline_preprocessed'];
data_file = dir([data_dir, data_name, '\*SZ*']);
funs_name = 'fdnoGRwrabrant_4D.nii';
save_dir = 'F:\fAPCS\workpipeline\SZImagingGenetic\results\functions\';

%% define ICs of RSNs
IC_dir = 'F:\fAPCS\workpipeline\SZImagingGenetic\3_preprocessing\_fMRI_ICA\';
IC25_name = 'rfMRI_ICA_d25.nii';
[mask_25,~,~,header] = y_ReadAll(fullfile(IC_dir, 'ukb_ICs', IC25_name));
IC25_components = 'rfMRI_GoodComponents_d25_v1.txt';
idx_m25 = load(fullfile(IC_dir, 'ukb_ICs', IC25_components));
mask_25_3mm = [];
for m = 1:length(idx_m25)
    dat = mask_25(:,:,:,idx_m25(m));
    NiftiWrite(dat, header, fullfile(IC_dir, 'ukb_ICs_2mm', ['rfMRI_ICA_d25_',int2str(idx_m25(m)), '_2mm.nii']));
    y_Reslice(fullfile(IC_dir, 'ukb_ICs_2mm', ['rfMRI_ICA_d25_',int2str(idx_m25(m)), '_2mm.nii']),fullfile(IC_dir, 'ukb_ICs_3mm', ['rfMRI_ICA_d25_',int2str(idx_m25(m)), '_3mm.nii']), [3 3 3], 1, 'ImageItself'); % 0: Nearest Neighbour. 1: Trilinear. 
    [temp,~,~,h] = y_ReadAll(fullfile(IC_dir, 'ukb_ICs_3mm', ['rfMRI_ICA_d25_',int2str(idx_m25(m)), '_3mm.nii'])); 
    mask_25_3mm = cat(4, mask_25_3mm, temp);
end
NiftiWrite(mask_25_3mm, h, fullfile(IC_dir, 'ukb_ICs', 'rfMRI_ICA_d25_3mm.nii'));

IC100_name = 'rfMRI_ICA_d100.nii';
[mask_100,~,~,header] = y_ReadAll(fullfile(IC_dir, 'ukb_ICs', IC100_name));
IC100_components = 'rfMRI_GoodComponents_d100_v1.txt';
idx_m100 = load(fullfile(IC_dir, 'ukb_ICs', IC100_components));
mask_100_3mm = [];
for m = 1:length(idx_m100)
    dat = mask_100(:,:,:,idx_m100(m));
    NiftiWrite(dat, header, fullfile(IC_dir, 'ukb_ICs_2mm', ['rfMRI_ICA_d100_',int2str(idx_m100(m)), '_2mm.nii']));
    y_Reslice(fullfile(IC_dir, 'ukb_ICs_2mm', ['rfMRI_ICA_d100_',int2str(idx_m100(m)), '_2mm.nii']),fullfile(IC_dir, 'ukb_ICs_3mm', ['rfMRI_ICA_d100_',int2str(idx_m100(m)), '_3mm.nii']), [3 3 3], 1, 'ImageItself'); % 0: Nearest Neighbour. 1: Trilinear. 
    [temp,~,~,h] = y_ReadAll(fullfile(IC_dir, 'ukb_ICs_3mm', ['rfMRI_ICA_d100_',int2str(idx_m100(m)), '_3mm.nii'])); 
    mask_100_3mm = cat(4, mask_100_3mm, temp);
end
NiftiWrite(mask_100_3mm, h, fullfile(IC_dir, 'ukb_ICs', ['rfMRI_ICA_d100_3mm.nii']));

%% read nifti data
mask_25 = spm_read_vols(spm_vol(fullfile(IC_dir, 'ukb_ICs', 'rfMRI_ICA_d25_3mm.nii')));
mask_25 = reshape(mask_25, [], size(mask_25, 4));
mask_100 = spm_read_vols(spm_vol(fullfile(IC_dir, 'ukb_ICs', 'rfMRI_ICA_d100_3mm.nii')));
mask_100 = reshape(mask_100, [], size(mask_100, 4));

sub_name = cell(length(data_file),1);
nodeamps25 = zeros(length(data_file), size(idx_m25,2));
net25 = zeros(length(data_file), size(idx_m25,2)*(size(idx_m25,2)-1)*0.5);
nodeamps100 = zeros(length(data_file), size(idx_m100,2));
net100 = zeros(length(data_file), size(idx_m100,2)*(size(idx_m100,2)-1)*0.5);

for sub = 1:length(data_file)
    sub_name{sub,1} = data_file(sub).name;
    gunzip(fullfile(data_dir, data_name, data_file(sub).name, [funs_name, '.gz']));
    data_4D = spm_read_vols(spm_vol(fullfile(data_dir, data_name, data_file(sub).name, funs_name)));
    data_2D = reshape(data_4D, [], size(data_4D,4));
    data_2D(isnan(data_2D)) = 0;
    delete(fullfile(data_dir, data_name, data_file(sub).name, funs_name));
    
    [tc_25, ~] = icatb_dual_regress(data_2D, mask_25); % y - Observations in columns (Voxels by time points), X - design matrix (Voxels by components)
    nodeamps25(sub, :) = std(tc_25);
    temp = corr(tc_25);
    temp = temp - diag(diag(temp));
    temp = 0.5 * log((1 + temp)./(1- temp));
    net25(sub, :) = squareform(temp);
    
    [tc_100, ~] = icatb_dual_regress(data_2D, mask_100);
    nodeamps100(sub, :) = std(tc_100);
    temp = corr(tc_100);
    temp = temp - diag(diag(temp));
    temp = 0.5 * log((1 + temp)./(1- temp));
    net100(sub, :) = squareform(temp);
end

%% save data
save_25 = [net25, nodeamps25];
save_100 = [net100, nodeamps100];

names_25 = cell(1, size(save_25, 2));
names_100 = cell(1, size(save_100, 2));
for i = 1:size(names_25, 2)
    if i <= size(nodeamps25,2)
        names_25{1,i} = ['NODEamps25 ', num2str(i)];
    else
        names_25{1,i} = ['NET25 ', num2str(i-size(nodeamps25,2))];
    end
end
for i = 1:size(names_100, 2)
    if i <= size(nodeamps100,2)
        names_100{1,i} = ['NODEamps100 ', num2str(i)];
    else
        names_100{1,i} = ['NET100 ', num2str(i-size(nodeamps100,2))];
    end
end
names_25 = ['ID', names_25];
save_25 = [sub_name, num2cell(save_25)];
save_25 = [names_25; save_25];
xlswrite(fullfile(save_dir, [data_idx, '_ICs_25.xlsx']), save_25);
names_100 = ['ID', names_100];
save_100 = [sub_name, num2cell(save_100)];
save_100 = [names_100; save_100];
xlswrite(fullfile(save_dir, [data_idx, '_ICs_100.xlsx']), save_100);
