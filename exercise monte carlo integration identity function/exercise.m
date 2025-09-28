% Exercise - Monte Carlo integration of identity function

%% 1. Aim of the exercise
% The aim of this exercise is to understand Monte Carlo integration by
% estimating the area under the identity function on the interval from 0 to
% 1 using random samples.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Define the number of samples 

% 3.1. Clear workspace and memory
clear;

% 3.2. Number of random samples to estimate the integral
N_samples = 1000;

%% 4. Monte Carlo integration: Estimating integrals via sampling

% 4.1. Generate random samples from Uniform[0,1]
uniform_samples = random('Uniform',0,1,[N_samples 1]);

% 4.2. MC estimate of the integral of f(x) = x over [0,1]
integral_estimate = mean(uniform_samples);

%% 5. The true value of the integral of the identity function over [0,1]

% 5.1. Set the true value of the integral
integral_true_value = 0.5; 

%% 6. Geometric interpretation of Monte Carlo integration

% 6.1 Create (x, y) points uniformly in the unit square
x_samples = uniform_samples; % x-coordinates sampled uniformly in [0,1]
y_samples = random('Uniform',0,1, ...
    [N_samples 1]); % y-coordinates sampled independently in [0,1]

% 6.2. Identify points below the diagonal y < x
below_diagonal = y_samples < x_samples; 

% 6.3. Visualizing area estimation via Monte Carlo sampling
figure
hold on
scatter(x_samples,y_samples, ...
    'DisplayName','All samples') % All samples
scatter(x_samples(below_diagonal),y_samples(below_diagonal), ...
    'DisplayName','Samples below diagonal'); % Samples below diagonal
plot([0 1],[0 1], ...
    'DisplayName','y = x') % Diagonal line y = x
title('Fig. 1. Monte Carlo estimation of area under f(x) = x')
xlabel('x')
ylabel('y')
legend('show')
hold off

%% 7. Convergence behavior of the Monte Carlo integral estimate

% 7.1. Track how the estimate evolves with more samples
convergence_integral_estimate = cumsum(uniform_samples)./(1:N_samples)';

% 7.2. Track how the mean squared error (MSE) evolves with more samples
convergence_MSE = cumsum((uniform_samples-integral_true_value).^2) ...
    ./(1:N_samples)';

% 7.3. Theoretical benchmark
theoretical_error_decay = 1./sqrt(1:N_samples);

%% 8. Visualize convergence of the Monte Carlo estimate and its error

% 8.1. Plot how the Monte Carlo estimate converges to the true value
figure
hold on
plot(1:N_samples,convergence_integral_estimate, ...
    'DisplayName','MC integral estimate');
yline(integral_true_value, ...
    'DisplayName','True integral value');
title('Fig. 2. Convergence of Monte Carlo estimate');
xlabel('Number of samples (draws)');
ylabel('Integral estimate');
legend('show');
hold off

% 8.3. Plot how the MSE decreases with more samples (log-log scale)
figure
hold on
loglog(1:N_samples,convergence_MSE, ...
    'DisplayName','MSE of MC estimate');
loglog(1:N_samples,theoretical_error_decay.^2, ...
    'DisplayName','Theoretical MSE decay');
title('Fig. 3. Log-log convergence of MC estimation error')
xlabel('Number of samples (log scale)')
ylabel('Mean Squared Error (log scale)')
legend('show')
hold off
