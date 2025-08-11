% Exercise - Monte Carlo integration of identity function

%% 1. Aim of the exercise
% To understand Monte Carlo (MC) integration by estimating the area under
% the identity function on the interval from 0 to 1 using random samples.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Define the number of samples 

% 3.1. Clear workspace and memory
clear;

% 3.2. Set number of random samples for estimation
N_samples = 1000; % Number of random draws from the target distribution

%% 4. Monte Carlo integration: Estimating integrals via sampling

% 4.1. Generate random samples from a uniform distribution on [0, 1]
uniform_samples = random('Uniform',0,1,[N_samples,1]);

% 4.2. MC estimate of the integral of f(x) = x over [0,1]
integral_estimate = mean(uniform_samples);

% Averaging f(x) = x over random values in [0,1] gives an estimate of the
% integral.

%% 5. The true value of the integral of the identity function over [0,1]

integral_true_value = 0.5; 

%% 6. Geometric interpretation of MC integration

% We sample random points (x,y) uniformly in the unit square [0,1] × [0,1].
% Each x-value (from x_samples) is treated as a fixed input to the function
% f(x) = x. Each y-value (from y_samples) is an independent random vertical
% coordinate. We do not compute f(x) directly. Instead, we test whether y <
% f(x). Points satisfying y < x lie below the curve and contribute to the
% area estimate. The fraction of such points approximates the integral of
% f(x) = x over [0,1].

% 6.1 Create (x, y) points uniformly in the unit square
x_samples = uniform_samples; % x-coordinates sampled uniformly in [0,1]
y_samples = random('Uniform',0,1, ...
    [N_samples,1]); % y-coordinates sampled independently in [0,1]

% 6.2. Identify points below the diagonal y < x
below_diagonal = y_samples < x_samples; 

% 6.3. Visualizing area estimation via Monte Carlo sampling
figure
hold on
scatter(x_samples,y_samples,'ro') % All samples
scatter(x_samples(below_diagonal),y_samples(below_diagonal), ...
    'o','MarkerEdgeColor','r','MarkerFaceColor', ...
    [0.9 0.9 0.9]); % Samples below diagonal
plot([0 1],[0 1]) % Diagonal line y = x
title('Fig. 1. Monte Carlo integration: Area under f(x) = x')
xlabel('x')
ylabel('y')
legend('All samples','Samples below diagonal','y = x')
hold off

%% 7. Convergence behavior of the MC integral estimate

% 7.1. Track how the estimate evolves with more samples
convergence_integral_estimate = cumsum(uniform_samples)./(1:N_samples)';

% cumsum computes the cumulative sum: [x1, x1+x2, x1+x2+x3, ...], where
% each xi is a sample from the uniform distribution. Dividing by [1, 2, 3,
% ...] gives the running mean at each step. This mimics adding one sample
% at a time and tracking how the estimate evolves. The running mean shows
% how our MC estimate of the integral evolves as more samples are added.
% Ideally, it should stabilize near the true value (0.5 for f(x) = x over
% [0,1]) due to the law of large numbers.

% 7.2. Track how the mean squared error (MSE) evolves with more samples
convergence_MSE = cumsum((uniform_samples-integral_true_value).^2) ...
    ./(1:N_samples)';

% Each term (xi-true_value)^2 measures the squared error of a sample
% relative to the true integral. cumsum accumulates these squared errors,
% and dividing by [1, 2, 3, ...] gives the running mean squared error. This
% shows how the average squared error of the MC estimate decreases as more
% samples are added. Ideally, the MSE should decay toward zero as sample
% size increases, reflecting improved accuracy.

% 7.3. Theoretical benchmark

theoretical_error_decay = 1./sqrt(1:N_samples);

% The MC error typically decreases proportionally to 1/sqrt(N) as the
% number of samples N increases. This line illustrates how fast we expect
% the error to decay with sample size. The reason: each sample has some
% variance, and averaging N samples reduces that variance by a factor of
% 1/N. The standard deviation (i.e., typical error) is the square root of
% the variance, hence it decays as 1/sqrt(N). This behavior arises from the
% Central Limit Theorem, which states that the distribution of the sample
% mean approaches a normal distribution as N increases, with variance
% decreasing as 1/N.

%% 8. Visualize convergence of the MC estimate and its error

% 8.1. Plot how the Monte Carlo estimate converges to the true value
figure
hold on
plot(1:N_samples,convergence_integral_estimate,'b','DisplayName', ...
    'MC integral estimate');
yline(integral_true_value,'DisplayName','True integral value');
title('Fig. 2. Convergence of Monte Carlo estimate');
xlabel('Number of samples (draws)');
ylabel('Integral estimate');
legend('show');
hold off

% 8.2. Plot how the MSE decreases with more samples
figure
hold on
plot(1:N_samples,convergence_MSE,'r','DisplayName','MSE of MC estimate');
plot(1:N_samples,theoretical_error_decay,'g','DisplayName', ...
    'Theoretical error decay');
title('Fig. 3. Convergence of MC estimation error metrics')
xlabel('Number of samples (draws)');
ylabel('Integral estimate error metrics');
legend('show');
hold off
