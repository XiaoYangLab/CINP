%% sort panss reduction after 6-week acute treatment of antipsychotic drugs
% Author: Yang Xiao, PKU, 2024
% @: xiaoyang9604@gmail.com

close all; 
clear all;
clc;

%% set working place
data_idx = ' '; 
data_dir = ' ';
data_name = ' ';
save_dir = ' ';

%% set factorization
pos_3 = 1:7;
neg_3 = 8:14;
gener = 15:30;

%% read panss data
[panss_dat, panss_txt] = xlsread(fullfile(data_dir, data_name),data_idx);
age = panss_dat(:,1);
sex = panss_dat(:,2);
baseline = panss_dat(:,3:32);
followup = panss_dat(:,33:62);

% 3 factor  
base_sum = sum(baseline,2);
base_pos_3 = sum(baseline(:,pos_3),2);
base_neg_3 = sum(baseline(:,neg_3),2);
base_gener = sum(baseline(:,gener),2);

follow_sum = sum(followup,2);
follow_pos_3 = sum(followup(:,pos_3),2);
follow_neg_3 = sum(followup(:,neg_3),2);
follow_gener = sum(followup(:,gener),2);

delta_sum = ((follow_sum-base_sum)./(base_sum-size(baseline,2))) .* 100;
delta_pos_3 = ((follow_pos_3-base_pos_3)./(base_pos_3-length(pos_3))) .* 100;
delta_neg_3 = ((follow_neg_3-base_neg_3)./(base_neg_3-length(neg_3))) .* 100;
delta_gener = ((follow_gener-base_gener)./(base_gener-length(gener))) .* 100;

%% save data
save_3 = [age, sex, delta_sum, delta_pos_3, delta_neg_3, delta_gener];
names_3 = [{'age'},{'sex'},{'allscores'},{'postive'},{'negative'},{'general'}];

names_3 = ['ID', names_3];
save_3 = [panss_txt(2:end,1), num2cell(save_3)];
save_3 = [names_3; save_3];

xlswrite(fullfile(save_dir, [data_idx, '_delta_panss_3.xlsx']), save_3);
