%   VIDEO: Permutation testing and shuffling

% number of 'trials' in each condition
n1 = 50;
n2 = 70; % note the trial inbalance!

% create data
data1 = randn(n1,1);
data2 = randn(n2,1) + .3; % note the mean offset!

% step 1: pool the data
alldata = cat(1, data1, data2);

% corresponding labels
truelabels = cat(1, ones(n1,1), 2*ones(n2,1) );

% compute the observed condition difference
idx_1 = (truelabels(:,1) == 1);
idx_2 = (truelabels(:,1) == 2);
true_conddif = mean(alldata(idx_2)) - mean(alldata(idx_1));


%%


% step 2: shuffle the data once
shuflabels = truelabels(randperm(length(truelabels)),1);


% step 3: compute the mean difference of the shuffled labels
idx_1 = (shuflabels(:,1) == 1);
idx_2 = (shuflabels(:,1) == 2);
shufconddif =   mean(alldata(idx_2)) - mean(alldata(idx_1));



%% creating a null-hypothesis (H0) distribution

% number of iterations for permutation testing
nIterations = 1000;

permvals = zeros(nIterations,1);


for permi=1:nIterations
    
    % steps 2 and 3 from above
    shuflabels = truelabels(randperm(length(truelabels)),1);
    idx_1 = (shuflabels(:,1) == 1);
    idx_2 = (shuflabels(:,1) == 2);
    permvals(permi) =  mean(alldata(idx_2)) - mean(alldata(idx_1));
end

% show the distribution
figure(1), clf, hold on
histogram(permvals,40)
plot([1 1]*true_conddif,get(gca,'ylim'),'r--','linew',3)
legend({'Shuffled';'True'})
set(gca,'xlim',[-1 1])
xlabel('Mean value')
ylabel('Count')


%% two methods for computing a p-value

% method p_z: compute the normalized distance of the observed p-value from
%             the distribution of H0 p-values

permmean = mean(permvals);
permstd  = std(permvals);

% formula for z-score
zdist = (true_conddif- permmean) / permstd;

% can convert to p-value
pval = 1-  normcdf(zdist);

title([ 'z_{dist}: ' num2str(zdist) ', p_z value: ' num2str(pval) ])


% method p_c
p_c = sum( permvals>true_conddif ) / nIterations;
title({[ 'z_{dist}: ' num2str(zdist) ', p_z value: ' num2str(pval) ]...
       [ 'p_c value: ' num2str(p_c) ]})
