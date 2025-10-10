% Exercise - Monte Carlo integration of a Gaussian function

%% 1. Aim of the exercise
% This exercise introduces Monte Carlo integration by estimating the area
% under a Gaussian curve using randomly sampled points. The goal is to 
% illustrate probabilistic approximation of definite integrals.

%% 2. Theory
% Refer to the accompanying PDF for theoretical background.

%% 3. Define the number of samples

% 3.1. Clear workspace and memory
clear;

% 3.2. Number of random samples to estimate the integral
N_samples = 1000; 

%% 4. Monte Carlo integration: Estimating the integral via sampling

% 4.1. Generate random samples from a standard normal distribution
standard_normal_samples = random('Normal',0,1,[N_samples 1]);

% 4.2. Evaluate the unnormalized Gaussian function at each sample
unnormalized_gaussian = exp(-0.5*standard_normal_samples.^2);

% 4.3. MC estimate of the integral 
integral_estimate = sqrt(2*pi)*mean(unnormalized_gaussian);

%% 5. Set the true value of the integral of the identity function

% Set the true value of the integral
integral_true_value = sqrt(pi); 

%% 6. Monte Carlo estimation error: Mean Squared Error (MSE)

% Compute the squared error of the Monte Carlo estimate
MSE = (integral_estimate-integral_true_value)^2;

%% 7. Convergence behavior of the MC integral estimate

% 7.1. Track how the estimate evolves with more samples
convergence_integral_estimate = ...
    cumsum(sqrt(2*pi).*unnormalized_gaussian)./(1:N_samples)';

% 7.2. Track how the mean squared error (MSE) evolves with more samples
% convergence_MSE = ...
%    cumsum((sqrt(2*pi).*unnormalized_gaussian- ...
%    integral_true_value).^2)./(1:N_samples)';
convergence_MSE = ...
    (convergence_integral_estimate-integral_true_value).^2;

%% 8. Convergence behavior of the Monte Carlo integral estimate

% 8.1. Theoretical benchmark
theoretical_error_decay = 1./sqrt(1:N_samples);

%% 9. Visualize convergence of the MC estimate and its error

% 9.1. Plot how the Monte Carlo estimate converges to the true value
figure
hold on
plot(1:N_samples,convergence_integral_estimate, ...
    'DisplayName','MC integral estimate');
yline(integral_true_value, ...
    'DisplayName','True integral value');
title('Fig. 1. Convergence of Monte Carlo estimate');
xlabel('Number of samples (draws)');
ylabel('Integral estimate');
legend('show');
hold off

% 9.2. Plot how the MSE decreases with more samples (log-log scale)
figure
hold on
loglog(1:N_samples,convergence_MSE, ...
    'DisplayName','MSE of MC estimate');
loglog(1:N_samples,theoretical_error_decay.^2, ...
    'DisplayName','Theoretical MSE decay');
title('Fig. 2. Log-log convergence of MC estimation error');
xlabel('Number of samples (log scale)');
ylabel('Mean Squared Error (log scale)');
legend('show');
hold off
