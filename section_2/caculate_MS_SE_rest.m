%% Calculate ms and se
%%
clc;clear;close all;

file_folder = dir('C:\Users\WANG\Desktop\rest_295\LR\REST_ROI_signals');        

file_folder = file_folder(3:end);            
num_brain = 100;                                     

for i = 1:length(file_folder)
    
    file_folder(i).name;
    load(['C:\Users\WANG\Desktop\rest_295\LR\REST_ROI_signals\',file_folder(i).name], 'REST_Net100_ROI');   
    
    
    for j = 1:num_brain    
        rest_Net100_ZROI = zscore(REST_Net100_ROI(:,j));      
        rest_all_LR(j,:,i) = rest_Net100_ZROI';                  
    end
   
end

for i = 1:length(file_folder)
    
    file_folder(i).name;
    load(['C:\Users\WANG\Desktop\rest_295\RL\REST_ROI_signals\',file_folder(i).name], 'REST_Net100_ROI');    % 下载静息态的bold信号
    
    
    for j = 1:num_brain    
        rest_Net100_ZROI = zscore(REST_Net100_ROI(:,j));     
        rest_all_RL(j,:,i) = rest_Net100_ZROI';                 
    end
    
end


rest_all(:,:,1:801) = rest_all_LR;
rest_all(:,:,802:1602) = rest_all_RL;

%% calculate the MS,SE and Kuramoto order paramater 

sub = 1602;
for i = 1:sub
    i
  
    M = int2str(i);
    
    
    ROI_246_RS = zscore(rest_all(:,:,i)');    
    
    signals_zs = ROI_246_RS;
   
    kop = xlz_kop(signals_zs');   
    
    
    bins_num = 30;
    [MS(i), SS(i), CS(i), SE(i), sample_failed] =  x_kop2sta(kop, bins_num);
    
    
    kop_sta_BN246(i).subj_ID = ROI_246_RS;   
    kop_sta_BN246(i).kop = kop;              
    kop_sta_BN246(i).mean_kop = MS(i);       
    kop_sta_BN246(i).min_kop = min(kop);     
    kop_sta_BN246(i).max_kop = max(kop);     
    kop_sta_BN246(i).std_kop = SS(i);        
    kop_sta_BN246(i).cv_kop = CS(i);         
    kop_sta_BN246(i).bins_num = bins_num;    
    kop_sta_BN246(i).entropy_kop = SE(i);    
    kop_sta_BN246(i).sample_failed = sample_failed;      
    
end
save(['rest_kuramoto_order_paramter.mat'],'kop_sta_BN246');





