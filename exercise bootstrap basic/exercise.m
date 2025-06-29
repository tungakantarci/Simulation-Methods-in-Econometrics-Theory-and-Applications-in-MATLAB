% Exercise - Understanding the method of basic bootsrap

%% 1. Aim of the exericse 
% To learn about the method of basic bootstrap.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Set values for the parameters of the simulation

% 3.1. Clear the memory 
clear;

% 3.2. Set the number of simulations as the number of bootstrap samples
N_sim = 1000;

%% 4. Create population data and draw a sample

% 4.1. Set the population size
N_obs_population = 10000;

% 4.2. Create the population data
data_population = random('Normal',4,5,[N_obs_population,1]);

% 4.3. Set the sample size
N_obs_sample = 500;

% 4.4. Draw a sample from the population
data_sample = datasample(data_population,N_obs_sample,'Replace',false);

%% 5. Draw (bootsrap) samples from the sample 

% 5.1. Preallocate vector to store the means of bootstrap samples
means_data_samples_bootstrap = NaN(N_sim,1);

% 5.2. Draw samples with replacement and compute the sample mean each time
for i = 1:N_sim
    data_samples_bootstrap = datasample(data_sample,N_obs_sample, ...
        'Replace',true);
    means_data_samples_bootstrap(i) = mean(data_samples_bootstrap);
end

%% 6. Draw samples from the population

% 6.1. Preallocate vector to store means of samples from the population
means_data_samples_population = NaN(N_sim,1);

% 6.2. Draw samples and compute the sample mean each time
for i = 1:N_sim
    data_samples_popuplation = datasample(data_population,N_obs_sample, ...
        'Replace',false);
    means_data_samples_population(i) = mean(data_samples_popuplation);
end

%% 7. Plot the PDFs of sample means from bootsrap and population sampling

% 7.1. Estimate the PDFs
[estimated_function_values_bootstrap,evaluation_points_bootstrap] = ...
    ksdensity(means_data_samples_bootstrap,'function','pdf');

[estimated_function_values_population,evaluation_points_population] = ...
    ksdensity(means_data_samples_population,'function','pdf');

% 7.2. Plot the PDFs
figure;
hold on
plot(evaluation_points_population,estimated_function_values_population);
plot(evaluation_points_bootstrap,estimated_function_values_bootstrap);
xlabel('Sample mean');
ylabel('Density');
title(['Fig. 1. PDFs of sample means based on bootstrap and ' ...
    'population sampling']);
legend('Population sampling','Bootstrap sampling');
hold off

%% 8. Plot the CDFs of sample means from bootsrap and population sampling

% 8.1. Estimate the CDFs
[estimated_function_values_bootstrap,evaluation_points_bootstrap] = ... 
    ksdensity(means_data_samples_bootstrap,'function','cdf');

[estimated_function_values_population,evaluation_points_population] = ...
    ksdensity(means_data_samples_population,'function','cdf');

% 8.2. Plot the CDFs
figure;
hold on
plot(evaluation_points_population,estimated_function_values_population);
plot(evaluation_points_bootstrap,estimated_function_values_bootstrap);
title(['Fig. 2. CDFs of sample means based on bootsrap and ' ...
    'population sampling']);
xlabel('Sample mean');
ylabel('Cumulative distribution');
legend('Population sampling','Bootstrap sampling');
hold off
