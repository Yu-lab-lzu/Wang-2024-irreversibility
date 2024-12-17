%% MS of different groups

clear;close all;clc;
load('ms_se_rest.mat')
mean_x = mean(x)

load('sub_middle_sup_260.mat','sub_260_1')
mean_M_ms_1 = mean(x(sub_260_1));

load('sub_middle_sup_260.mat','sub_260_2')
mean_M_ms_2 = mean(x(sub_260_2));

load('sub_middle_sup_260.mat','sub_260_3')
mean_M_ms_3 = mean(x(sub_260_3));

load('sub_middle_sup_260.mat','sub_260_4')
mean_M_ms_4 = mean(x(sub_260_4));

load('sub_middle_sup_260.mat','sub_260_5')
mean_M_ms_5 = mean(x(sub_260_5));

load('sub_middle_sup_260.mat','sub_260_6')
mean_M_ms_6 = mean(x(sub_260_6));

load('sub_middle_sup_260.mat','sub_260_7')
mean_M_ms_7 = mean(x(sub_260_7));

load('sub_middle_sup_260.mat','sub_260_8')
mean_M_ms_8 = mean(x(sub_260_8));

load('sub_middle_sup_260.mat','sub_260_9')
mean_M_ms_9 = mean(x(sub_260_9));

load('sub_middle_sup_260.mat','sub_260_10')
mean_M_ms_10 = mean(x(sub_260_10));

load('sub_middle_sup_260.mat','sub_260_11')
mean_M_ms_11 = mean(x(sub_260_11));

load('sub_middle_sup_260.mat','sub_260_12')
mean_M_ms_12 = mean(x(sub_260_12));

a = [mean_M_ms_1 mean_M_ms_2 mean_M_ms_3 mean_M_ms_4  mean_M_ms_5   mean_M_ms_6   ...
    mean_M_ms_7  mean_M_ms_8  mean_M_ms_9  mean_M_ms_10  mean_M_ms_11 mean_M_ms_12];

