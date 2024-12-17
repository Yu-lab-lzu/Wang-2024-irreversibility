%% Define state transition

clear;close all;clc;

load('Data_processed_motor','cell_IDs','spike_times','stimuli_lengths','task_3','num_brain')


num_cellSample = 100;
num_cells = num_brain;
num_stimuli = 1;
sub = 270;


dt = 4000;

n = 6;

cell_samples = cell(1,num_cellSample);
transitions = cell(num_stimuli,num_cellSample);
num_repeats = [sub];


for i = 1:num_cellSample
    i
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

        % 记录内容
        transitions{j,i} = transitions_temp;      % 每个刺激重复的每一次都放在矩阵 transitions 中
    end
    toc
end

save('transitions_motor_N6_4000.mat','n', 'num_stimuli','num_cellSample',...
    'num_repeats', 'stimuli_lengths', 'cell_samples', 'transitions','cell_IDs','spike_times');






