% Exercise - Monte Carlo integration of logarithmic function 

%% 1. Aim of the exercise
% This script estimates the number of prime numbers less than or equal to n
% using Monte Carlo integration to approximate the logarithmic integral,
% which serves as a smooth estimate of the prime-counting function π(n).

%% 2. Theory
% Refer to the accompanying PDF for theoretical background.

%% 3. Simulation Setup

% 3.1. Clear workspace and memory
clear;

% 3.2. Set number of random samples for Monte Carlo estimation
N_samples = 1000;  % Number of random draws from Uniform(2, n)

% 3.3. Define evaluation points for prime counting and integration
N_values = 100:100:10000;

% Create a vector of n values ranging from 100 to 10,000 in increments of
% 100. For each n, we will estimate the logarithmic integral using Monte
% Carlo integration, and compute the actual number of primes smaller than
% or equal to n This allows us to compare the approximation against true
% prime counts across a broad range.

%% 4. Monte Carlo Integration: Estimating the Logarithmic Integral

% 4.1. Preallocate vector for Monte Carlo estimates
MC_log_integral = NaN(size(N_values));  

% 4.2. Preallocate vector for actual prime counts
prime_counts = NaN(size(N_values));     

% 4.3. Loop over each n value and perform Monte Carlo estimation
for i = 1:length(N_values)
    
    % Current value of n
    N = N_values(i);
    
    % Draw random samples uniformly from [2, N]
    uniform_samples = random('Uniform',2,N,[N_samples,1]);
    
    % Apply transformation: f(x) = 1/log(x)
    log_transformation = 1./log(uniform_samples);
    
    % Estimate integral using sample mean
    mean_value = mean(log_transformation);
    MC_log_integral(i) = (N-2)*mean_value;
    
    % Compute actual number of primes up to N
    prime_counts(i) = length(primes(N));
end

%% 5. Visualization: Compare Integral Estimates with True Prime Counts

hold on
plot(N_values, MC_log_integral,'o','DisplayName','MC estimate', ...
    'MarkerFaceColor','red');
plot(N_values, prime_counts,'o','DisplayName','Actual prime count', ...
    'MarkerFaceColor','blue');
title(['Fig. 1. Monte Carlo estimation of logarithmic integral for ' ...
    'prime counting']);
xlabel('n');
ylabel('Number of primes smaller than or equal to n');
legend('show')
hold off
