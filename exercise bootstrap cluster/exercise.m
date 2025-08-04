% Exercise - Understanding the method of cluster bootstrap

%% 1. Aim of the exercise  
% To learn how to estimate regression coefficients and their standard
% errors using two complementary methods: bootstrap resampling of clusters
% and analytical cluster-robust estimation using the sandwich formula.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Load data

% 3.1. Clear the memory 
clear;

% 3.2. Load dataset from CSV file
data = readtable('data.csv');

%% 4. Set the number of simulations as the number of bootstrap samples

% 4.1. Set the number of bootstrap samples
N_sim = 1000;

% 4.2. Define the number of coefficient estimates to be simulated
N_k = 4;

%% 5. Determine the number of clusters

% 5.1. Obtain cluster identifiers
cluster_identifiers = data.country;  

%. 5.2. Identify unique clusters
unique_clusters = unique(cluster_identifiers);

% 5.3. Count the number of clusters
N_clusters = numel(unique_clusters);

%% 6. Draw (bootsrap) samples from the original sample 

% 6.1. Preallocate vector to store coefficient estimates
B_hats_data_samples_boot = NaN(N_k,N_sim);

% 6.2 Bootstrap: Sample clusters with replacement and compute coef. estimates
for i = 1:N_sim
    % Randomly sample clusters with replacement
    idx = randi(N_clusters,N_clusters,1);
    sampled_clusters = unique_clusters(idx);

    % Initialize cell array to store data for the selected cluster
    cluster_chunks = cell(numel(sampled_clusters),1);

    % Extract rows corresponding to the sampled cluster
    for j = 1:numel(sampled_clusters)
        cluster_chunks{j} = data(strcmp(data.country,sampled_clusters(j)),:);
    end

    % Combine data from all selected clusters into one table
    resampled_data = vertcat(cluster_chunks{:});

    % Define response variable
    y = table2array(resampled_data(:,4));

    % Create intercept term
    x_0 = ones(size(resampled_data,1),1);

    % Create the design matrix
    X = [x_0,table2array(resampled_data(:,[3,5,6]))];

    % Run robust OLS regression and store coefficient estimates
    LSS = exercisefunctionlssrobust(y,X);
    B_hats_data_samples_boot(:,i) = LSS.B_hat;
end

%% 7. Estimate cluster-robust SEs using the bootsrap resampling

% 7.1. Mean of the bootstrap sampling distribution of coefficient estimates
bootstrap_mean = mean(B_hats_data_samples_boot,2);

% 7.2. Standard deviation of the bootstrap sampling distribution (bootstrap SEs)
bootstrap_SE = std(B_hats_data_samples_boot,0,2);

%% 8. Define variables for full-sample OLS regression

% 8.1. Define response variable
y = table2array(data(:,4));

% 8.2. Create intercept term
x_0 = ones(size(data,1),1);

% 8.3. Create the design matrix
X = [x_0,table2array(data(:,[3,5,6]))];

%% 9. Obtain OLS statistics
LSS = exercisefunctionlssrobust(y,X);

%% 10. Preallocate matrices for cluster-robust SE estimation

% 10.1. Initialize cluster-level meat matrix (X_c'*u_c*u_c'*X_c)
M_c = zeros(N_k,N_k);

% 10.2. Initialize sum of meat matrices across clusters
M_sum = zeros(N_k,N_k);

% 10.3. Initialize cluster-level design matrix cross-product (X_c'*X_c)
XTX_c = zeros(N_k,N_k);

% 10.4. Initialize sum of design matrix cross-products across clusters
XTX_sum = zeros(N_k,N_k);

%% 11. Compute cluster-level contributions to SE estimator

% 11.1. Loop over each cluster to compute meat and bread components
for i = 1:N_clusters
    % Identify current cluster
    cluster_name = unique_clusters(i);

    % Find rows corresponding to the current cluster
    cluster_rows = strcmp(data.country,cluster_name);

    % Extract residuals for the current cluster
    error_c = LSS.u_hat(cluster_rows);

    % Extract design matrix rows for the current cluster
    X_c = X(cluster_rows,:);

    % Compute cluster-level meat matrix
    M_c = X_c'*(error_c*error_c')*X_c;

    % Accumulate meat matrices
    M_sum = M_sum+M_c;

    % Compute cluster-level design matrix cross-product
    XTX_c = X_c'*X_c;

    % Accumulate design matrix cross-products
    XTX_sum = XTX_sum+XTX_c;
end

%% 12. Estimate cluster-robust SEs using the analytical sandwich formula

% 12.1. Define the sample size
N_obs = size(data,1);

% 12.2. Apply finite sample correction for cluster-robust SE estimator
cluster_df_correction = (N_obs-1)/(N_obs-N_k)*(N_clusters/(N_clusters-1));

% 12.3. Compute cluster-robust SEs using the sandwich formula
cluster_robust_SE = cluster_df_correction*sqrt(diag(inv(XTX_sum)*M_sum*inv(XTX_sum))); % Sandwich formula: sqrt(diag(correction_factor*(X'X)^(-1)*meat*(X'X)^(-1)))
