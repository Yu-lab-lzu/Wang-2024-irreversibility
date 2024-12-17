%% Define state transition
clear;close all;clc;

for gg = 1:1:12

    gg

    load(['Data_processed_rest_1200_',num2str(gg),'.mat'],'cell_IDs','spike_times','stimuli_lengths','task_3','num_brain')

    num_cellSample = 1000;                  
    num_cells = num_brain;                
    num_stimuli = 1;                       
    sub = 260;
    
    dt = 5000;

    n = 5;

    cell_samples = cell(1,num_cellSample);                
    transitions = cell(num_stimuli,num_cellSample);       
    num_repeats = [sub];

    for i = 1:num_cellSample           
        i
        gg
        tic
        % Pick random group of regions 
        cells = datasample(1:num_cells, n,'Replace',false);     

        %  Record group of regions
        cell_samples{i} = cells;       

        % Loop over different stimuli
        for j = 1:num_stimuli        

            % Record transitions:
            transitions_temp = cell(num_repeats(j),1);        

            % Loop over different stimuulus repeats   
            for k = 1:num_repeats(j)         

                if isempty(cell2mat(spike_times{j}(k)))
                    T = zeros(2^n, n+1);
                else
                    % Compute transitions:
                    T = transitions_slidingWindow_variableLengths(spike_times{j}(k),...
                        cell_IDs{j}(k), cells, dt, cell2mat(stimuli_lengths{j}(k)));
                end
                
                transitions_temp{k} = sparse(T);

            end

           transitions{j,i} = transitions_temp;     
        end
        toc
    end

    save(['transitions_rest_N5_5000_1200_',num2str(gg),'.mat'],'n', 'num_stimuli','num_cellSample',...
        'num_repeats', 'stimuli_lengths', 'cell_samples', 'transitions','cell_IDs','spike_times');


end






