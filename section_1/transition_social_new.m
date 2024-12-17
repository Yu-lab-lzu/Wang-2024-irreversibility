%% Define state transition


clear;close all;clc;

load('Data_processed_social.mat','cell_IDs','spike_times','stimuli_lengths','task_3','num_brain')


num_cellSample = 100;          
num_cells = num_brain;   
num_stimuli = 1;               
sub = 270;

dt = 4000;   

%
n = 6;

cell_samples = cell(1,num_cellSample);                
transitions = cell(num_stimuli,num_cellSample);       
num_repeats = [sub];


for i = 1:num_cellSample           

    i
    tic
    
    cells = datasample(1:num_cells, n,'Replace',false);     
    
    
    cell_samples{i} = cells;       
    
    
    for j = 1:num_stimuli         
        
        
        transitions_temp = cell(num_repeats(j),1);        
        
        
        for k = 1:num_repeats(j)         
            
            if isempty(cell2mat(spike_times{j}(k)))
                T = zeros(2^n, n+1);
            else
            
            T = transitions_slidingWindow_variableLengths(spike_times{j}(k),...
                cell_IDs{j}(k), cells, dt, cell2mat(stimuli_lengths{j}(k)));
            end
            

            transitions_temp{k} = sparse(T);
            
        end
        
        
        transitions{j,i} = transitions_temp;     
    end
    toc
end


save('transitions_social_N6_4000.mat','n', 'num_stimuli','num_cellSample',...
    'num_repeats', 'stimuli_lengths', 'cell_samples', 'transitions','cell_IDs','spike_times');









