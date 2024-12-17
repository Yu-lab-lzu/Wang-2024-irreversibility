%% compute irreversibility
clear;close all;clc;
for gg = 1:1:12

    % Load state transitions:
    data_struct = load(['transitions_rest_N5_5000_1200_',num2str(gg),'.mat']);

    % Transitions for each repeat of each states and each region group:
    transitions = data_struct.transitions;

    % Number of repeats for each stimulus:
    num_repeats = data_struct.num_repeats;

    % Lengths of different repeats of different state (in ms):
    stimuli_lengths = data_struct.stimuli_lengths;

    % Number of state and samples:
    num_stimuli = data_struct.num_stimuli;
    num_samples = data_struct.num_cellSample;

    % Size of region groups:
    n = data_struct.n;

    %% 

    % Data fractions to consider:
    fracs = [0.9:-0.1:0.5];
    num_fracs = length(fracs);

    % Number of data subsamples for each data fraction:
    num_dataSamples = 100;

    % Initial step size for irreversibility algorithm:
    step_size_init = .5;

    % Number of iterations of irreversibility algorithm:
    num_steps = 300;

    % Transpose indices:
    inds_trans = transpose_indices_multipartite(n);

    % Observables for different constraints:
    Os = cell(n-1,1);
    for i = 1:(n-1)
        Os{i} = binary_constraints_multipartite(n, i);
    end

    % List of region samples:
    cell_samples = data_struct.cell_samples;

    % Record number of transitions:
    num_trans = zeros(num_stimuli, num_samples);

    % Record irreversibilities and inifinite-data estimates:
    irreversibilities = zeros(num_stimuli, num_samples, num_fracs + 1, num_dataSamples);
    irreversibilities_inf = zeros(num_stimuli, num_samples, num_dataSamples);

    % Record minimum irreversibilities of order k from 1 to N-1:
    irreversibilities_min = zeros(num_stimuli, n-1, num_samples, num_fracs + 1, num_dataSamples);
    irreversibilities_min_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);

    % Record proportions of minimum and interaction irreversibilities:
    props_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);
    props_int_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);

    % Record changes in state probabilities:
    change_stateProbs = zeros(num_stimuli, num_samples, 2^n, num_fracs + 1, num_dataSamples);
    change_stateProbs_inf = zeros(num_stimuli, num_samples, 2^n, num_dataSamples);

    % Loop over different state:
    for i = 1:num_stimuli   

        % Loop over different region samples:
        for j = 1:num_samples     

            tic
            gg
            % Copmute full transition matrix:
            T = zeros(2^n, n+1);
            T_shuffle = zeros(2^n, n+1);

            for k = 1:num_repeats(i)
                T = T + full(transitions{i,j}{k});
            end

            % Number of transitions:
            num_trans(i,j) = sum(T(:));

            % List all transitions:
            trans = [];

            for k = 1:2^n
                for l = 1:(n+1)
                    trans = [trans; repmat([k, l], T(k,l), 1)];
                end
            end

            % Compute irreversibility with full data (using pseudocount method):
            P = (T(:,1:n) + 1)/sum(T(:) + 1);
            irreversibilities(i,j,1,:) = sum(P(:).*log2(P(:)./P(inds_trans)));

            % Compute independent irreversibility with full data:
            p_1 = reshape(sum(Os{1}.*repmat(P, 1, 1, 2*n), [1 2]), 2, n);
            irreversibilities_min(i,1,j,1,:) = sum(p_1.*log2(p_1./flipud(p_1)), [1 2]);

            % Compute change in state probabilities:
            % P_full = convert_transitions_multipartite((T + 1)/sum(T(:) + 1));
            % change_stateProbs(i,j,:,1,:) = repmat(reshape(sum(P_full,1) - sum(P_full,2)', 1, 1, 2^n), 1, 1, 1, 1, num_dataSamples);

            % Keep track of transitions:
            Ps = zeros(2^n, n, num_fracs, num_dataSamples);

            % Loop over data subsamples:
            for k = 1:num_dataSamples

                % Initialize list of transitions:
                trans_temp = trans;

                % Loop over fractions of data:
                for l = 1:num_fracs

                    % Sample transitions:
                    inds = randsample(size(trans_temp,1), round(fracs(l)*num_trans(i,j)), false);
                    trans_temp = trans_temp(inds,:);
                    [trans_unique, ~, ic] = unique(trans_temp, 'rows');
                    T_temp = zeros(2^n, n + 1);
                    for m = 1:size(trans_unique, 1)
                        T_temp(trans_unique(m,1), trans_unique(m,2)) = sum(ic == m);
                    end

                    % Compute irreversibility (using pseudocounts):
                    P_temp = (T_temp(:,1:n) + 1)/sum(T_temp(:) + 1);
                    irreversibilities(i,j,l+1,k) = sum(P_temp(:).*log2(P_temp(:)./P_temp(inds_trans)));

                    Ps(:,:,l,k) = P_temp;

                    % Compute independent irreversibility:
                    p_1 = reshape(sum(Os{1}.*repmat(P_temp, 1, 1, 2*n), [1 2]), 2, n);
                    irreversibilities_min(i,1,j,l+1,k) = sum(p_1.*log2(p_1./flipud(p_1)), [1 2]);

                    % Compute change in state probabilities:
                    % P_full = convert_transitions_multipartite((T_temp + 1)/sum(T_temp(:) + 1));
                    % change_stateProbs(i,j,:,l+1,k) = reshape(sum(P_full,1) - sum(P_full,2)', 1, 1, 2^n);

                end

                % Extrapolate irreversibility to infinite data with linear fit:
                fit = polyfit([1, fracs].^(-1), reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
                irreversibilities_inf(i,j,k) = fit(2);

                % Extrapolate independent irreversibility to infinite data:
                fit_ind = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,1,j,:,k), 1, num_fracs + 1), 1);
                irreversibilities_min_inf(i,1,j,k) = fit_ind(2);
                fit_prop_ind = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,1,j,:,k), 1, num_fracs + 1)./...
                    reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
                props_inf(i,1,j,k) = fit_prop_ind(2);

              
            end

            % Print some things:
            i
            j
            toc

        end
    end
    clear Os transitions data_struct irreversibilities_min change_stateProbs_inf fit_prop_ind props_inf
    clear change_stateProbs trans_temp num_trans T trans irreversibilities_min_inf T_shuffle props_int_inf
    clear dd k j num_samples num_repeats num_dataSamples trans_unique P cell_samples
    save(['entropy_5000_N5_260/sub_entropy_rest_N5_5000_',num2str(gg),'.mat'])
    datetime
end






