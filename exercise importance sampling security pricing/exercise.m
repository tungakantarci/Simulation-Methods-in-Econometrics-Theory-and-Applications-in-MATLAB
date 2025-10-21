% Exercise - Understanding importance sampling with security pricing

%% 1. Aim of the exercise
% The aim of this exercise is to understand importance sampling in the 
% context of reducing variance of the integral estimate in a security 
% pricing example.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Clear the memory

% Clear the memory 
clear;

%% 4. Define dividend function

% 4.1. Define dividend function as anonymous function
d_fun = @(x) sin((2*pi)/21.*x).*cos((2*pi)/exp(1).*x);

% 4.2. Create sequence for plotting
x = -5:0.1:25;

% 4.3. Evaluate dividend function over sequence
d = d_fun(x);

%% 5. Visualize the shape of the dividend function

% Create figure and plot dividend function over x
figure
plot(x,d)
title('Fig. 1. Dividend function')
xlabel('x')
ylabel('d')
ylim([-2 2])

%% 6. Create the trajectory of the dividend random variable

% 6.1. Draw samples from Gamma distribution
X = random('Gamma',1.2,0.5,[1 101]);

% 6.2. Evaluate dividend realizations
d_X = abs(d_fun(X));

%% 7. Visualize the trajectory of the dividend random variable

% Create figure and plot dividend realizations
figure
plot(-5:0.1:5,d_X)
title('Fig. 2. Trajectory of the dividend random variable');
xlabel('Sample index');
ylabel('Dividend (currency units)');

%% 8. Compute the weighted importance sampling estimate

% 8.1. Set proposal distribution: Exponential with rate 2 (mean = 0.5)
proposal_mean = 0.5;

% 8.2. Define weight function w(x) = 2*x^0.2
w_fun = @(x) 2.*x.^0.2;

% 8.3. Define scaled dividend function g(x) = (1/R)*|d(x)|
g_fun = @(x) (1/1.05).*abs(d_fun(x));

% 8.4. Create vector of sample sizes from 1,000 to 100,000
N_sizes = 1000:1000:100000;

% 8.5. Preallocate vector for importance sampling estimates
weighted_is_estimate = zeros(size(N_sizes));

% 8.6. Preallocate vector for sample variances of estimates
sample_variance_expo = zeros(size(N_sizes));

% 8.7. Loop over each sample size to compute IS estimate and variance
for i = 1:length(N_sizes)
    
    % 8.7.1. Set current sample size
    N = N_sizes(i);
    
    % 8.7.2. Draw N samples from Exponential(mean = 0.5)
    Xp = random('Exponential',proposal_mean,[N 1]);
    
    % 8.7.3. Evaluate scaled dividend function g(X_i)
    g_vals = g_fun(Xp);
    
    % 8.7.4. Compute importance weights w(X_i)
    w_vals = w_fun(Xp);
    
    % 8.7.5. Normalize weights to sum to 1
    w_norm = w_vals/sum(w_vals);
    
    % 8.7.6. Compute weighted importance sampling estimate
    weighted_is_estimate(i) = sum(w_norm.*g_vals);
    
    % 8.7.7. Compute sample variance of the IS estimate
    sample_variance_expo(i) = (N/(N-1))* ...
        sum((w_norm.^2).*(g_vals-weighted_is_estimate(i)).^2);
end

%% 9. Visualize convergence behavior of importance sampling estimate

% Create figure and plot IS estimates across sample sizes
figure
plot(N_sizes,weighted_is_estimate,'DisplayName','IS estimate')
title('Fig. 3. Convergence of the IS estimate')
xlabel('Sample size')
ylabel('Estimated price (currency units)')
legend('show')

%% 10. Visualize SD behavior of importance sampling estimate

% Create figure for log-transformed standard deviation comparison
figure
log_sd = log(sqrt(sample_variance_expo)); % Compute log of sample SD from IS estimates
log_n_half = log(1./sqrt(N_sizes)); % Compute log of theoretical O(n^{-1/2}) convergence rate
log_n_inv = log(1./N_sizes); % Compute log of theoretical O(n^{-1}) convergence rate
hold on
plot(N_sizes,log_sd,'DisplayName','SD of IS estimate');
plot(N_sizes,log_n_half,'DisplayName','O(n^{-1/2})');
plot(N_sizes,log_n_inv,'DisplayName','O(n^{-1})');
title('Fig. 4. Standard deviation of importance sampling estimate');
xlabel('Sample size');
ylabel('log(Standard deviation)');
legend('show');
hold off
