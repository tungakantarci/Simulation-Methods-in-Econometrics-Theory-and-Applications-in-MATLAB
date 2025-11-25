% Exercise - Understanding importance sampling using tail probability

%% 1. Aim of the exercise

% The aim of this exercise is to undersand importance sampling as a
% technique for reducing the variance of Monte Carlo integration
% estimators.

%% 2. Theory

% Refer to the accompanying PDF file for the theory.

%% 3. Clear the memory

% Clear the memory 
clear;

%% 4. Define and plot tail probability region

% 4.1. Define evaluation grid centered at z_thresh
z_base = -5:0.1:5;

% 4.2. Define threshold for tail probability
z_thresh = 1.96;

% 4.3. Ensure z_thresh is included and grid is sorted
z = unique([z_base,z_thresh]); 

% 4.4. Define Standard Normal distribution
pd_base_norm = makedist('Normal','mu',0,'sigma',1);

% 4.5. Evaluate PDF over the grid
PDF_target_norm = pdf(pd_base_norm,z);

% 4.6. Region to highlight where z is at least z_thresh
z_tail = z >= z_thresh;

% 4.7. Plot Standard Normal PDF and highlight tail region
figure
hold on
plot(z,PDF_target_norm, ...
    'Color',[0.0000 0.0000 0.0000], ...
    'DisplayName','N(0,1)');
fill( ...
    [z(z_tail),fliplr(z(z_tail))], ...
    [PDF_target_norm(z_tail),zeros(size(PDF_target_norm(z_tail)))], ...
    [0.9000 0.9000 0.9000], ...
    'EdgeColor','none', ...
    'DisplayName','Tail probability: P(Z >= z\_thresh)');
xline(z_thresh,':', ...
    'DisplayName','z threshold');
title('Fig. 1. Standard Normal tail probability beyond z\_thresh');
xlabel('z');
ylabel('PDF');
legend('show');
hold off

%% 5. Compute true tail probability

% 5.1. Define mean of the Normal distribution
mu = 0; 

% 5.2. Define standard deviation of the Normal distribution
sigma = 1;

% 5.3. Evaluate tail probability beyond z_thresh as true integral value
int_true_value = 1-cdf('Normal',z_thresh,mu,sigma);

%% 6. Define number of samples

% Number of random samples to estimate the integral
N_samples = 100000;

%% 7. Monte Carlo integration using standard sampling

% 7.1. Draw samples from the Standard Normal distribution
Z_samples = random('Normal',mu,sigma,[N_samples 1]);

% 7.2. Estimate tail probability using indicator function
int_est_mc = mean(Z_samples >= z_thresh);

% 7.3. Compute MSE of the Monte Carlo estimator
MSE_mc = 1/(N_samples*(N_samples-1)) * ...
    sum(((Z_samples >= z_thresh)-int_est_mc).^2);

%% 8. Importance sampling integration using Normal proposal

% 8.1. Define location parameter of Normal proposal
mu = z_thresh;

% 8.2. Define scale parameter of Normal proposal
sigma = 0.5;

% 8.3. Draw samples from Normal proposal
norm_samples = random('Normal',mu,sigma,[N_samples 1]);

% 8.4. Compute importance weights: target density/proposal density
weights_norm = ...
    normpdf(norm_samples,0,1)./ ... 
    normpdf(norm_samples,mu,sigma); 

% 8.5. Estimate tail probability via importance sampling
int_est_is_norm = mean(weights_norm.*(norm_samples >= z_thresh));

% 8.6. Compute MSE of the importance sampling estimator
MSE_norm = var(weights_norm.*(norm_samples >= z_thresh))/N_samples;

%% 9. Importance sampling integration using Truncated Normal proposal

% 9.1. Define location parameter of Normal proposal
mu = z_thresh;

% 9.2. Define scale parameter of Normal proposal
sigma = 0.5;

% 9.3. Define base Normal distribution for truncation
pd_base_trunc_norm = makedist('Normal','mu',mu,'sigma',sigma);

% 9.4. Define Trunc. Normal proposal: N(z_thresh,0.5) on [z_thresh,inf)
pd_prop_trunc_norm = truncate(pd_base_trunc_norm,z_thresh,inf);

% 9.5. Draw samples from the Truncated Normal distribution
trunc_norm_samples = random(pd_prop_trunc_norm,[N_samples 1]);

% 9.6. Compute importance weights
weights_trunc_norm = ...
    normpdf(trunc_norm_samples,0,1) ./ ...
    (normpdf(trunc_norm_samples,mu,sigma) ./ ...
    (1 - normcdf(z_thresh,mu,sigma)));

% 9.7. Estimate tail probability via importance sampling
int_est_is_trunc_norm = mean(weights_trunc_norm); 

% 9.8. Compute MSE of the importance sampling estimator
MSE_trunc_norm = var(weights_trunc_norm)/N_samples;

%% 10. Importance sampling integration using Gamma proposal

% 10.1. Define shape parameter of Gamma proposal
alpha = 5;

% 10.2. Define scale parameter of Gamma proposal
theta = 0.5;

% 10.3. Draw samples from Gamma proposal
gamma_samples = random('Gamma',alpha,theta,[N_samples 1]);

% 10.4. Select samples in the integration domain: z >= z_thresh
gamma_selected = gamma_samples(gamma_samples >= z_thresh);

% 10.5. Compute importance weights
weights_gamma = ...
    normpdf(gamma_selected,0,1)./ ... % Standard Normal
    gampdf(gamma_selected,alpha,theta); % Gamma

% 10.6. Estimate tail probability via importance sampling
int_est_is_gamma = mean(weights_gamma);

% 10.7. Compute MSE of the importance sampling estimator
MSE_gamma = var(weights_gamma)/length(weights_gamma);

%% 11. Compare target and proposal PDFs

% 11.1. Define full evaluation grid for symmetric distributions
x = -5:0.1:10;

% 11.2. Define positive evaluation grid for Gamma support
x_pos = 0:0.1:10;

% 11.3. Define Truncated evaluation grid for Truncated Normal
x_trunc = z_thresh:0.1:10;

% 11.4. Create Standard Normal distribution object
pd_target_norm = makedist('Normal','mu',0,'sigma',1);

% 11.5. Evaluate Standard Normal PDF over full grid
PDF_target_norm = pdf(pd_target_norm,x);

% 11.6. Create Normal proposal distribution centered at z_thresh
pd_prop_norm = makedist('Normal','mu',z_thresh,'sigma',0.5);

% 11.7. Evaluate Normal proposal PDF over full grid
PDF_prop_norm = pdf(pd_prop_norm,x);

% 11.8. Create base Normal distribution for truncation
pd_base_trunc_norm = makedist('Normal','mu',z_thresh,'sigma',0.5);

% 11.9. Truncate base Normal distribution to [z_thresh,inf)
pd_prop_trunc_norm = truncate(pd_base_trunc_norm,z_thresh,inf);

% 11.10. Evaluate Truncated Normal PDF over truncated grid
PDF_prop_trunc_norm = pdf(pd_prop_trunc_norm,x_trunc);

% 11.11. Create Gamma proposal distribution
pd_prop_gamma = makedist('Gamma',5,0.5);

% 11.12. Evaluate Gamma proposal PDF over positive grid
PDF_prop_gamma = pdf(pd_prop_gamma,x_pos);

% 11.13. Create comparison plot
figure
hold on
plot(x,PDF_target_norm, ...
    'Color',[0.0000 0.0000 0.0000], ...
    'DisplayName','Target: N(0,1)');
xline(z_thresh,':', ...
    'Color',[0.0000 0.0000 0.0000], ...
    'DisplayName','z\_thresh');
plot(x,PDF_prop_norm, ...
    'Color',[0.8500 0.3250 0.0980], ...
    'DisplayName','Proposal: N(z\_thresh,0.5)');
plot(x_trunc,PDF_prop_trunc_norm, ...
    'Color',[0.9290 0.6940 0.1250], ...
    'DisplayName','Proposal: TruncNorm(z\_thresh,0.5)');
plot(x_pos,PDF_prop_gamma, ...
    'Color',[0.4940 0.1840 0.5560], ...
    'DisplayName','Proposal: Gamma(5,0.5)');
title('Fig. 2. Target vs. Proposal PDFs');
xlabel('z');
ylabel('PDF');
legend('show');
hold off

%% 12. Compare likelihood ratios

% 12.1. Define evaluation grid
x = -5:0.1:10;

% 12.2. Define evaluation grid for truncated proposals
x_trunc = z_thresh:0.1:10;

% 12.3. Define evaluation grid for positive-support proposals
x_pos = 0:0.1:10;

% 12.4. Create Standard Normal distribution object
pd_target_norm = makedist('Normal','mu',0,'sigma',1);

% 12.5. Evaluate Standard Normal PDF over full grid
PDF_target_norm = pdf(pd_target_norm,x);

% 12.6. Evaluate Standard Normal PDF over truncated grid
PDF_target_trunc_norm = pdf(pd_target_norm,x_trunc);

% 12.7. Evaluate Standard Normal PDF over positive grid
PDF_target_pos_norm = pdf(pd_target_norm,x_pos);

% 12.8. Create Normal proposal distribution
pd_prop_norm = makedist('Normal','mu',z_thresh,'sigma',0.5);

% 12.9. Evaluate Normal proposal PDF over full grid
PDF_prop_norm = pdf(pd_prop_norm,x);

% 12.10. Create base Normal distribution for truncation
pd_base_trunc_norm = makedist('Normal','mu',z_thresh,'sigma',0.5);

% 12.11. Truncate base Normal distribution to [z_thresh, inf)
pd_prop_trunc_norm = truncate(pd_base_trunc_norm,z_thresh,inf);

% 12.12. Evaluate Truncated Normal proposal PDF over truncated grid
PDF_prop_trunc_norm = pdf(pd_prop_trunc_norm,x_trunc);

% 12.13. Create Gamma proposal distribution
pd_prop_gamma = makedist('Gamma',5,0.5);

% 12.14. Evaluate Gamma proposal PDF over positive grid
PDF_prop_gamma = pdf(pd_prop_gamma,x_pos);

% 12.15. Create likelihood ratio plot
figure
hold on
plot(x,PDF_target_norm./PDF_prop_norm, ...
    'Color',[0.8500 0.3250 0.0980], ...
    'DisplayName','Norm(z\_thresh,0.5)');
plot(x_trunc,PDF_target_trunc_norm./PDF_prop_trunc_norm, ...
    'Color',[0.9290 0.6940 0.1250], ...
    'DisplayName','TruncNorm(z\_thresh,0.5)');
plot(x_pos,PDF_target_pos_norm./PDF_prop_gamma, ...
    'Color',[0.4940 0.1840 0.5560], ...
    'DisplayName','Gamma(5,0.5)');
xline(z_thresh,':', ...
    'Color',[0.0000 0.0000 0.0000], ...
    'DisplayName','z\_thresh');
xlim([0 5]); % xlim([-5 10]) is full domain, including symmetric region
ylim([0 1]); % Rescale y-axis to reveal structure
title('Fig. 3. Likelihood ratios: Target / Proposal');
xlabel('z');
ylabel('Likelihood ratio');
legend('show');
hold off

%% 13. Convergence behavior of importance sampling estimators

% Here we conduct an experiment to evaluate the convergence behavior of
% different estimators by comparing their performance at two representative
% sample sizes: 1,000 and 100,000. Rather than embedding the full
% estimation procedure in a for loop, fixed samples are used to assess how
% each estimator behaves under varying conditions, focusing on qualitative
% consistency across runs.
