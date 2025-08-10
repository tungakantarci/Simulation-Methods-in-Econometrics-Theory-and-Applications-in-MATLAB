% Exercise - Understanding the method of block bootstrap

%% 1. Aim of the exericse 
% To understand how to correctly estimate the sampling distribution of a
% statistic (e.g., the sample mean) when data are serially correlated, such
% as in time-series settings. Standard bootstrap methods assume independent
% observations and can lead to biased inference. We use the block bootstrap
% method to account for dependence across observations and compare its
% performance to sampling directly from the population.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Set values for the parameters of the simulation

% 3.1. Clear the memory
clear;

% 3.2. Set the number of simulations as the number of bootstrap samples
N_sim = 5000;

%% 4. Generate population data

% 4.1. Set the population size
N_obs_pop = 100000;

% 4.2. Set AR(1) persistence parameter
phi = 0.7;

% 4.3. Generate white noise innovations
varepsilon = random('Normal',0,1,[N_obs_pop,1]);

% 4.4. Simulate AR(1) time‑series data
data_pop = filter(1,[1-phi],varepsilon);

%% 5. Plot AR(1) time series
figure;
plot(1:1000,data_pop(1:1000));
xlabel('t');
ylabel('X_t');
title('Fig. 1. First 1000 Points of the AR(1) Series');

%% 6. Draw samples from the population

% 6.1. Set the sample size
N_obs_sample = 5000;

% 6.2. Preallocate matrix to store samples
data_samples_pop = NaN(N_obs_sample,N_sim);

% 6.3. Preallocate vector to store sample means
means_data_samples_pop = NaN(N_sim,1);

% 6.4. Draw samples from the population and compute the sample mean each time
for i = 1:N_sim
    start_idx = random('Discrete Uniform',N_obs_pop-N_obs_sample+1,[1,1]); 
    sample_i = data_pop(start_idx:start_idx+N_obs_sample-1);
    data_samples_pop(:,i) = sample_i;
    means_data_samples_pop(i) = mean(sample_i);
end

%% 7. Generate sample data

% 7.1. Randomly select a starting index for sampling
start_idx = random('Discrete Uniform',N_obs_pop-N_obs_sample+1,[1,1]); % Ensures sample fits within bounds of data_pop

% 7.2. Draw one contiguous sample from the population
data_sample = data_samples_pop(:,1); % Use the first sample drawn from the population

% 7.3. Compute optimal block length
block_length = ceil(N_obs_sample^(1/5)); % Rule-of-thumb for block bootstrap

% 7.4. Compute how many full blocks can be formed
N_blocks = floor(N_obs_sample/block_length); % Maximum number of full blocks

% 7.5. Trim sample to fit an exact number of full blocks
trimmed_sample = data_sample(1:N_blocks*block_length);  % Ensures reshape works cleanly

% 7.6. Reshape trimmed sample into non-overlapping blocks
blocks = reshape(trimmed_sample,block_length,N_blocks)';

%% 8. Draw (bootsrap) samples from the initial sample 

% 8.1. Preallocate vector to store (bootsrap) sample means
means_data_samples_boot = NaN(N_sim,1);

% 8.2. Draw samples of k blocks from the original sample and compute the sample mean each time
for i = 1:N_sim
    data_samples_boot = datasample(blocks,N_blocks,1,'Replace',true);
    data_samples_boot_flat = data_samples_boot(:);
    means_data_samples_boot(i) = mean(data_samples_boot_flat);
end

%% 9. Plot the estimated PDFs
figure;
hold on
ksdensity(means_data_samples_boot,'function','pdf');
ksdensity(means_data_samples_pop,'function','pdf');
title(['Fig. 2. PDFs of sample means based on bootstrap and ' ...
    'population sampling']);
xlabel('Sample mean');
ylabel('Density');
legend('Bootstrap sample means','Population sample means');
hold off

%% 10. Plot the estimated CDFs
figure;
hold on
ksdensity(means_data_samples_boot,'function','cdf');
ksdensity(means_data_samples_pop,'function','cdf');
title(['Fig. 3. CDFs of sample means based on bootsrap and ' ...
    'population sampling']);
xlabel('Sample mean');
ylabel('Cumulative distribution');
legend('Bootstrap sample means','Population sample means');
hold off
