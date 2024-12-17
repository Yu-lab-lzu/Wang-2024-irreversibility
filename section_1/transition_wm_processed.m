%% The data processing stage before calculating the state transition matrix T of the wm task
 
clear;close all;clc;

%%
file_folder = dir('C:\Users\WANG\Desktop\ENTROPY\wm\network_100_glm');        
file_folder = file_folder(3:end);         
num_brain = 100;

for i = 1:length(file_folder)
    
    file_folder(i).name;
    load(['C:\Users\WANG\Desktop\ENTROPY\wm\network_100_glm\',file_folder(i).name,'\WM_Net_100_ROI.mat'], 'WM_Net_100_ROI');

    for j = 1:num_brain
        wm_Net100_ZROI = zscore(WM_Net_100_ROI(:,j));   
        wm_all(j,:,i) = wm_Net100_ZROI';
    end

end
wm_all(:,177:end,:) = [];

%%

sub = 270;       
m = 1;           
L = 176;        

task_3 = cell(1,m);     
task_3{1,1} = wm_all;

for i = 1:sub     
    for j = 1:num_brain    
        
        wm_all_pca = wm_all(j,:,i);           
        
        std_bold(j,i) = std(wm_all_pca(:));   
        
    end
end


stimuli_time = cell(1,m);       
stimuli_IDs = cell(1,m);        
for u = 1:m                     
    
    stimuli_time{1,u} = cell(1,sub);
    stimuli_IDs{1,u} = cell(1,sub);
    
end


for u = 1:m        
    for i = 1:sub      
        for j = 1:num_brain  
            
            [a1 a2] = find(task_3{1,u}(j,:,i) > std_bold(j,i)); 
  
            bold = task_3{1,u}(j,:,i);   
            
            b = cell(0,0);
            num = [];    
            ct = 1;        
            head = 1;      
            tail = 1;      
            
            while (ct < numel(a2))               
                
                head = ct;   
                ct = ct+1;   
                
                while (ct <= numel(a2) && (a2(ct) - a2(ct-1))==1)  
                    ct = ct+1;
                end
                tail = ct - 1;          
                b = [b;a2(head:1:tail)'];      
                num = [num;tail-head+1];   
            end
            
            if(tail<numel(a2))                          
                b = [b;a2(tail+1)'];                    
                num = [num; 1];                        
            end
            %
            if (numel(a2) == 1 )    
                b = [b,a2(1)];
            end
            
            if isempty(a2)
                
                stimuli_time{1,u}{j,i} = [];        
                stimuli_IDs{1,u}{j,i} = [];   
                
            else
                for k = 1:length(b)     
                    
                    A = b{k,1};
                    
                    if A(1) == 1;     
                        AA(k) = 0;
                    else
                        AA(k) = A(1)-((task_3{1,u}(j,A(1),i) - std_bold(j,i))./ (task_3{1,u}(j,A(1),i) - task_3{1,u}(j,A(1)-1,i)));
                    end
                    
                end
                
                AA(find(AA == 0)) = [];
                stimuli_time{1,u}{j,i} = AA;        
                clear AA
                
            end
        end
    end
end

for i = 1:m                       
    for j = 1:sub              
        for k = 1:num_brain       
            
        num_spike_point{1,i}{k,j} = numel(stimuli_time{1,i}{k,j});
            
        end
        num_spike_point{1,i}{k+1,j} = sum(cell2mat(num_spike_point{1,i}(:,j)));
    end
end



cell_IDs = cell(1,m);
spike_times = cell(1,m);
stimuli_lengths = cell(1,m);             
for i = 1:m
    
    cell_IDs{1,i} = cell(1,sub);       
    spike_times{1,i} = cell(1,sub);
    stimuli_lengths{1,i} = cell(1,sub);
    A1{1,i} = cell(1,sub);
    
end


%% 
t = 720;
A = [];
B = [];

for i = 1:m         
    
    tic
    
    for j = 1:sub       
        A = [];
        B = [];
        
        for k = 1:num_brain    
            
            if isempty(cell2mat(stimuli_time{1,i}(k,j)))
                k = k+1;         
            else
                A1{1,i}(1,j) = mat2cell([A,cell2mat(stimuli_time{1,i}(k,j))],1);  
                A = cell2mat(A1{1,i}(1,j));
            end
            
        end
        
        
        [a,b] = sort(cell2mat(A1{1,i}(1,j)));     
        
        for k = 1:num_brain
            
            B = [B,k*ones(1,cell2mat(num_spike_point{1,i}(k,j)))];            
            
        end


        C = B(b);                                     
        cell_IDs{1,i}{1,j} =  C;
        spike_times{1,i}{1,j} = t*roundn(a,-2);

        stimuli_lengths{1,i}{1,j} = L.*t;          
        
    end
    toc
    
end

clear a1 a2 num bold ct head tail A a b B C

save('Data_processed_WM','cell_IDs','spike_times','stimuli_lengths','task_3','num_brain')



