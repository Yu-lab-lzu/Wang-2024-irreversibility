%%
clear;close all;clc;

data_struct = load('transitions_rest_N6_4000');


transitions = data_struct.transitions;


num_repeats = data_struct.num_repeats;

stimuli_lengths = data_struct.stimuli_lengths;

num_stimuli = data_struct.num_stimuli;
num_samples = data_struct.num_cellSample;

n = data_struct.n;

%% 

fracs = [0.98:-0.04:0.74];
num_fracs = length(fracs);


num_dataSamples = 100;

step_size_init = .5;


num_steps = 300;


inds_trans = transpose_indices_multipartite(n);


Os = cell(n-1,1);
for i = 1:(n-1)
    Os{i} = binary_constraints_multipartite(n, i);
end


cell_samples = data_struct.cell_samples;


num_trans = zeros(num_stimuli, num_samples);

irreversibilities = zeros(num_stimuli, num_samples, num_fracs + 1, num_dataSamples);
irreversibilities_inf = zeros(num_stimuli, num_samples, num_dataSamples);

irreversibilities_min = zeros(num_stimuli, n-1, num_samples, num_fracs + 1, num_dataSamples);
irreversibilities_min_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);

props_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);
props_int_inf = zeros(num_stimuli, n-1, num_samples, num_dataSamples);

change_stateProbs = zeros(num_stimuli, num_samples, 2^n, num_fracs + 1, num_dataSamples);
change_stateProbs_inf = zeros(num_stimuli, num_samples, 2^n, num_dataSamples);


for i = 1:num_stimuli    
    
    
    for j = 1:num_samples     
        
        tic
        
        
        T = zeros(2^n, n+1);
        T_shuffle = zeros(2^n, n+1);
        
        for k = 1:num_repeats(i)
            T = T + full(transitions{i,j}{k});
        end
        
        num_trans(i,j) = sum(T(:));
        
        trans = [];
        
        for k = 1:2^n
            for l = 1:(n+1)
                trans = [trans; repmat([k, l], T(k,l), 1)];
            end
        end
        
        P = (T(:,1:n) + 1)/sum(T(:) + 1);
        irreversibilities(i,j,1,:) = sum(P(:).*log2(P(:)./P(inds_trans)));
        
        p_1 = reshape(sum(Os{1}.*repmat(P, 1, 1, 2*n), [1 2]), 2, n);
        irreversibilities_min(i,1,j,1,:) = sum(p_1.*log2(p_1./flipud(p_1)), [1 2]);
        
        P_full = convert_transitions_multipartite((T + 1)/sum(T(:) + 1));
        change_stateProbs(i,j,:,1,:) = repmat(reshape(sum(P_full,1) - sum(P_full,2)', 1, 1, 2^n), 1, 1, 1, 1, num_dataSamples);
        
        Ps = zeros(2^n, n, num_fracs, num_dataSamples);
        
        for k = 1:num_dataSamples
        
            trans_temp = trans;
            
            for l = 1:num_fracs
                
                inds = randsample(size(trans_temp,1), round(fracs(l)*num_trans(i,j)), false);
                trans_temp = trans_temp(inds,:);
                [trans_unique, ~, ic] = unique(trans_temp, 'rows');
                T_temp = zeros(2^n, n + 1);
                for m = 1:size(trans_unique, 1)
                    T_temp(trans_unique(m,1), trans_unique(m,2)) = sum(ic == m);
                end
                
                P_temp = (T_temp(:,1:n) + 1)/sum(T_temp(:) + 1);
                irreversibilities(i,j,l+1,k) = sum(P_temp(:).*log2(P_temp(:)./P_temp(inds_trans)));
                
                Ps(:,:,l,k) = P_temp;
                
                p_1 = reshape(sum(Os{1}.*repmat(P_temp, 1, 1, 2*n), [1 2]), 2, n);
                irreversibilities_min(i,1,j,l+1,k) = sum(p_1.*log2(p_1./flipud(p_1)), [1 2]);
                
                P_full = convert_transitions_multipartite((T_temp + 1)/sum(T_temp(:) + 1));
                change_stateProbs(i,j,:,l+1,k) = reshape(sum(P_full,1) - sum(P_full,2)', 1, 1, 2^n);
                  
            end
            
            fit = polyfit([1, fracs].^(-1), reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
            irreversibilities_inf(i,j,k) = fit(2);
            
            fit_ind = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,1,j,:,k), 1, num_fracs + 1), 1);
            irreversibilities_min_inf(i,1,j,k) = fit_ind(2);
            fit_prop_ind = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,1,j,:,k), 1, num_fracs + 1)./...
                reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
            props_inf(i,1,j,k) = fit_prop_ind(2);
            
            for l = 1:(2^n)
                fit_dP = polyfit([1, fracs].^(-1), reshape(change_stateProbs(i,j,l,:,k), 1, num_fracs + 1), 1);
                change_stateProbs_inf(i,j,l,k) = fit_dP(2);
            end
            
        end
 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        if mean(irreversibilities_inf(i,j,:)) > 2*std(irreversibilities_inf(i,j,:))
            
            P_temp = P;
            
            for k = (n-1):-1:2
                
                [P_temp, S_temp, ~] = min_irreversibility_multipartite(Os{k}, P_temp, num_steps, step_size_init);
                irreversibilities_min(i,k,j,1,:) = S_temp;
                
            end
            
            for k = 1:num_dataSamples
                
                for l = 1:num_fracs
                    
                    P_temp = Ps(:,:,l,k);
                    
                    for m = (n-1):-1:2
                        
                        [P_temp, S_temp, ~] = min_irreversibility_multipartite(Os{m}, P_temp, num_steps, step_size_init);
                        irreversibilities_min(i,m,j,l+1,k) = S_temp;
                        
                    end
                    
                end
                
                for l = 2:(n-1)
                    
                    fit_irreversibility = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,l,j,:,k), 1, num_fracs + 1), 1);
                    irreversibilities_min_inf(i,l,j,k) = fit_irreversibility(2);
                    
                    fit_prop = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,l,j,:,k), 1, num_fracs + 1)./...
                        reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
                    props_inf(i,l,j,k) = fit_prop(2);
                    
                    fit_prop_int = polyfit([1, fracs].^(-1), reshape(irreversibilities_min(i,l,j,:,k) - irreversibilities_min(i,l-1,j,:,k), 1, num_fracs + 1)./...
                        reshape(irreversibilities(i,j,:,k), 1, num_fracs + 1), 1);
                    props_int_inf(i,l,j,k) = fit_prop_int(2);
                    
                end
            end
        end
        
        
        i
        j
        toc
        
    end
end

 
props_int_inf(:,1,:,:) = props_inf(:,1,:,:);


sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
sig_min = (abs(mean(irreversibilities_min_inf, 4)) > 2*std(irreversibilities_min_inf, [], 4)).*sign(mean(irreversibilities_min_inf, 4));
sig_change_stateProbs = (abs(mean(change_stateProbs_inf, 4)) > 2*std(change_stateProbs_inf, [], 4)).*sign(mean(change_stateProbs_inf, 4));


inds_sig = cell(num_stimuli,1);
inds_sig_min = cell(num_stimuli, n-1);

for i = 1:num_stimuli
    
    inds_sig{i} = find(sig(i,:) > 0);
    
    for j = 1:(n-1)
        inds_sig_min{i,j} = find(sig_min(i,j,:) > 0);
    end
end

% Save results:
save('irreversibility_rest_N6_4000', 'num_stimuli', 'num_repeats',...
    'num_dataSamples', 'fracs', 'stimuli_lengths', 'cell_samples', 'num_trans',...
    'irreversibilities', 'irreversibilities_min', 'irreversibilities_inf',...
    'irreversibilities_min_inf', 'props_inf', 'props_int_inf', 'sig', 'sig_min',...
    'inds_sig', 'inds_sig_min');

