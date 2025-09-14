% Exercise - Understanding the CLT using simulation

%% 1. Aim of the exercise
% To learn the CLT using simulation.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Set the parameters of the simulation

% 3.1. Clear the memory
clear;

% 3.2. Define the population size
N_obs_population = 10000;

% 3.3. Define alternative sample sizes
N_obs_sample = [2,15,30,90];

% 3.4. Define the number of samples
N_samples = size(N_obs_sample,2);

% 3.5. Define the number of simulated samples 
N_sim = 1000;

%% 4. An exponential random variable
% Demonstrating the CLT using an exponential variable.

%% 4.1. Define the population

% 4.1.1. Define Lambda
Lambda = 1;

% 4.1.2. Define the population 
  population = random('Exponential',Lambda,[N_obs_population 1]);
% population = random('Uniform',0,2,[N_obs_population 1]);

%% 4.2 Plot the frequency distribution of the exponential variable

% 4.2.1. Create the plot
figure
histogram(population);
title('Fig. 1. Histogram of the exponential variable');
ylabel('Frequncy');
xlabel('Values');

%% 4.3. Plot the sampling distribution of the sample mean

% 4.3.1. Preallocate an array to store means of samples
means_samples = NaN(N_sim,N_samples);

% 4.3.2. Draw random samples from the population and take their mean
for i = 1:N_sim
    for j = 1:N_samples
        sample = randsample(population,N_obs_sample(j),true);
        means_samples(i,j) = mean(sample);
    end    
end

% 4.3.3. Create the plot
colors = ['b','r','g','m'];
figure
hold on
for j = 1:N_samples
    [estimated_function_values_j,evaluation_points_j] = ...
        ksdensity(means_samples(:,j));
    plot(evaluation_points_j,estimated_function_values_j, ...
        colors(mod(j-1,length(colors))+1));
    title('Fig. 2. Sampling distribution of the sample mean');
    ylabel('Probability density');
    xlabel('Values');
end
legend_labels = arrayfun(@(x) sprintf('Sample size = %d', ...
    N_obs_sample(x)),1:N_samples,'UniformOutput',false);
legend(legend_labels);
hold off

%% 4.4. Speed of convergence of the sampling distribution

% 4.4.1. Define the theoretical skewness of the normal distribution
theoretical_skewness = 0;

% 4.4.2. Preallocate matrix to store skewness values
means_samples_skewness = NaN(1,N_samples);

% 4.4.3. Calculate the skewness of the sampling distribution
for j = 1:N_samples
    means_samples_skewness(1,j) = skewness(means_samples(:,j));
end

% 4.4.4. Define the absolute difference
abs_dif = abs(means_samples_skewness-theoretical_skewness);

% 4.4.5. Create the plot
figure
hold on
for j = 1:N_samples
    scatter(N_obs_sample(j),abs_dif(j),1000,'Marker','.', ...
        'DisplayName',sprintf('Sample size = %d',N_obs_sample(j)));
end
ylim([0 1.5]);
title(['Fig. 3. The difference between the skewness of ' ...
    'the sampling and theoteical distributions']);
ylabel('Difference in skewness');
xlabel('Sample size');
legend('show');
hold off
