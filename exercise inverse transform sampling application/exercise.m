% Exercise – Understanding inverse transform sampling

%% 1. Aim of the exercise

% The aim of this exercise is to understand the method of inverse transform
% sampling using the exponential and truncated exponential distributions.

%% 2. Theory

% Refer to the accompanying PDF for the theory of inverse transform
% sampling.

%% 3. Clear workspace

% Clear
clear;

%% 4. Sampling setup

% 4.1. Sample size
N = 1000;

%% 5. Sampling from the exponential distribution
% 
% 5.1 Generate N realizations from Uniform(0,1)
U = random('Uniform',0,1,[N 1]); 

% 5.2. Rate parameter of the exponential distribution
lambda = 1;

% 5.3. Apply the inverse CDF of Exp(lambda)
ExponentialRealization = -(1/lambda)*log(1-U);

%% 6. Histogram of exponential sample

% Plot
figure
hold on
histogram(ExponentialRealization,50,'Normalization','pdf', ...
    'EdgeColor',[0.0 0.2 0.4], ...
    'FaceColor',[0.0 0.2 0.4]);
ylim([0.0 1.2])
xlim([-0.5 9])
title('Fig. 1. Exponential distribution: Inverse-transform sample');
hold off

%% 7. MATLAB's built-in exponential generator

% Generate exponential sample using MATLAB's built‑in function
X_builtin = random('Exponential',1/lambda,[N 1]);

%% 8: Histogram of MATLAB's built‑in exponential sample

% Plot
figure
hold on
histogram(X_builtin,50,'Normalization','pdf', ...
    'EdgeColor',[0.0 0.2 0.4], ...
    'FaceColor',[0.0 0.2 0.4]);
ylim([0.0 1.2])
xlim([-0.5 9])
title('Fig. 2. Exponential distribution: Built-in sample');
hold off

%% 9. Sampling from the truncated exponential

% 9.1. Truncation point
a = 2;

% 9.2. Draw random numbers between 0 and 1
U = random('Uniform',0,1,[N 1]);

% 9.3. Create the sample that begins at a and extends to infinity
X_trunc = a-log(1-U);

%% 10. Histogram of truncated exponential sample

% Plot
figure
hold on
histogram(X_trunc,50,'Normalization','pdf', ...
    'EdgeColor',[0.0 0.2 0.4], ...
    'FaceColor',[0.0 0.2 0.4]);
ylim([0.0 1.2])
xlim([-0.5 9])
title(['Fig. 3. Truncated exponential distribution: ' ...
    'Inverse-transform sample']);
hold off
