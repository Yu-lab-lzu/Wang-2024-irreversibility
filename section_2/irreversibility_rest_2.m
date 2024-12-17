%%  compute irreversibility

load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_1.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

    
end
rest_sub_entropy_1 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_2.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

   
end
rest_sub_entropy_2 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))



%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_3.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

    
end
rest_sub_entropy_3 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))



%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_4.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

    
end
rest_sub_entropy_4 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_5.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

    
end
rest_sub_entropy_5 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))

%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_6.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

end
rest_sub_entropy_6 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_7.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

    
end
rest_sub_entropy_7 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_8.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

    
end
rest_sub_entropy_8 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


%
load('entropy_5000_N5_260_2\sub_entropy_rest_N5_5000_9.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

    
end
rest_sub_entropy_9 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_10.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

    
end
rest_sub_entropy_10 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_11.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

   
end
rest_sub_entropy_11 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


%
load('entropy_5000_N5_260\sub_entropy_rest_N5_5000_12.mat')
sig = (abs(mean(irreversibilities_inf, 3)) > 2*std(irreversibilities_inf, [], 3)).*sign(mean(irreversibilities_inf, 3));
inds_sig = cell(num_stimuli,1);

for i = 1:num_stimuli

    inds_sig = find(sig(i,:) > 0);

end
rest_sub_entropy_12 = mean(mean(squeeze(irreversibilities_inf(i,inds_sig,:)),2))


b = [rest_sub_entropy_1 rest_sub_entropy_2 rest_sub_entropy_3 rest_sub_entropy_4  rest_sub_entropy_5,...
     rest_sub_entropy_6  rest_sub_entropy_7 ,...
    rest_sub_entropy_8  rest_sub_entropy_9  rest_sub_entropy_10,...
     rest_sub_entropy_11 rest_sub_entropy_12];

save('data_500_entropy_260.mat','a','b')

