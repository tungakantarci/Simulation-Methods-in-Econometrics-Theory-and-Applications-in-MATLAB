%% Exercise - Monte Carlo integration of a Gaussian function

%% 1. Aim of the exercise
% Use Monte Carlo (MC) integration to estimate the area under a Gaussian
% curve. Instead of computing the integral analytically, approximate it
% using random sampling.

%% 2. Theory
% Refer to the accompanying PDF for theoretical background.

%% 3. Define the number of samples

% 3.1. Clear workspace and memory
clear;

% 3.2. Set number of random samples for estimation
N_samples = 1000; % Number of random draws from the target distribution

%% 4. Monte Carlo integration: Estimating integrals via sampling

% 4.1. Generate random samples from a standard normal distribution
standard_normal_samples = random('Normal',0,1,[N_samples,1]);

% 4.2. Evaluate the unnormalized Gaussian function at each sample point
unnormalized_gaussian = exp(-0.5*standard_normal_samples.^2);

% This computes the height of the Gaussian curve at each sampled x-value.

% 4.3. Estimate the integral by averaging the function values and scaling
integral_estimate = sqrt(2*pi)*mean(unnormalized_gaussian);

%% 5. Convergence behavior of the MC Estimate

% 5.1. Track how the estimate evolves with increasing sample size
convergence_integral_estimate = ...
    cumsum(sqrt(2*pi).*unnormalized_gaussian)./(1:N_samples)';

% Running average of scaled function values shows convergence toward the
% true value.

% 5.2. Compute running mean squared error (MSE) of the estimate
convergence_MSE = ...
    cumsum((sqrt(2*pi).*unnormalized_gaussian-sqrt(pi)).^2) ./ ...
    (1:N_samples)';

% Computes the running MSE between the scaled Gaussian estimate and
% sqrt(pi). The numerator accumulates squared errors; the denominator
% normalizes by sample count.

%% 6. Visualize convergence of the MC estimate and error metrics

% 6.1. Plot how the Monte Carlo estimate converges to the true value
figure
hold on
plot(1:N_samples,convergence_integral_estimate,'DisplayName', ...
    'MC integral estimate');
yline(sqrt(pi),'DisplayName','True integral value');
title('Fig. 1. Convergence of Monte Carlo estimate');
xlabel('Number of samples (draws)');
ylabel('Estimate');
legend('show');
hold off

% 6.2. Plot how the MSE decreases with increasing sample size
figure
hold on
plot(1:N_samples,convergence_MSE,'DisplayName','MSE of MC estimate');
plot(1:N_samples,1./sqrt(1:N_samples),'DisplayName', ...
    'Theoretical error decay');
title('Fig. 2. Convergence of MSE of MC estimate');
xlabel('Number of samples (draws)');
ylabel('Mean Squared Error');
legend('show');
hold off
