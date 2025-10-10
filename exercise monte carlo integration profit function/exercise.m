% Exercise – Understanding Monte Carlo integration using a profit function

%% 1. Aim of the exercise
% To estimate the total profit over a simulated volume using Monte Carlo
% (MC) integration. The profit depends on time, seasonal effects,
% inflation, and input variables such as x and y.

%% 2. Theory
% Refer to the accompanying PDF file for the theoretical background on the
% profit model and MC integration.

%% 3. Define the number of samples

% 3.1. Clear workspace and memory
clear;

% 3.2. Set number of random samples
N_samples = 1000000;

%% 4. Monte Carlo sampling of input variables

% 4.1. Define total volume over time and space
volume_v = 365*100*100;

% 4.2. Generate random x values
x_sample = 100*random('Uniform',0,1,[N_samples 1]);

% 4.3. Generate continuous y values
y_sample = 100*random('Uniform',0,1,[N_samples 1]);

% 4.4. Generate continuous time values (days)
t_sample = 365*random('Uniform',0,1,[N_samples 1]);

%% 5. Compute profit components

% 5.1. Inflation adjustment: compound interest model
inflation = (1+0.02/365).^t_sample;

% 5.2. Price function: adjusted for inflation and input costs
p_function = 5.*inflation-(1/200).*x_sample-(1/300).*y_sample;

% 5.3. Cost function: increases with inputs and time
c_function = (2+0.015*x_sample+0.01*y_sample).*(1+t_sample./1000);

% 5.4. Seasonal effect: models periodic variation over time
seasonal_effect = (1/3)*cos((2*pi/365).*t_sample-pi/6)+1;

% 5.5. Growth factor: exponential increase over time
growth = exp(t_sample./1000);

% 5.6. Quantity function: depends on inputs and seasonal effects
q_function = (80-0.05.*x_sample.*seasonal_effect-0.08.*y_sample) ...
    .*growth;

% 5.7. Profit per sample: revenue minus cost, scaled by quantity
pi_function = (p_function-c_function).*q_function;

%% 6. Estimate total profit

% Monte Carlo estimate of total profit
profit_estimate = volume_v*mean(pi_function);

%% 7. Construct confidence interval for integral estimate

% 7.1. Compute sample variance of profit evaluations
sample_variance = var(pi_function);

% 7.2. Compute SE of the profit estimate
profit_estimate_SE = volume_v*sqrt(sample_variance/N_samples);

% 7.3. Construct 95% CI for the profit estimate
CI_lower = profit_estimate-1.96*profit_estimate_SE;
CI_upper = profit_estimate+1.96*profit_estimate_SE;

%% 8. Diagnostic: Theoretical MSE of the Monte Carlo integral

% 8.1. Estimate theoretical MSE of the Monte Carlo integral (scaled)
MSE_scaled = volume_v^2*sample_variance/N_samples;

% 8.2. Estimate normalized MSE (per unit volume)
MSE_normalized = sample_variance/N_samples;
