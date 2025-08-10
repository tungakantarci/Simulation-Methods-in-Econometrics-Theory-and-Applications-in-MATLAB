% Exercise - Understanding the method of paired bootstrap

%% 1. Aim of the exercise
% To understand how to correctly estimate standard errors of regression
% coefficients when the error variance is not constant (i.e.,
% heteroskedasticity is present). Standard inference methods can be
% misleading in such cases. We use the paired bootstrap method, which
% resamples observations as (y, X) pairs, to account for heteroskedasticity
% and compare its performance to sampling directly from the population.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Set values for the parameters of the simulation

% 3.1. Clear the memory
clear;

% 3.2. Set the number of simulations as the number of bootstrap samples
N_sim = 5000;

%% 4. Generate population data

% 4.1. Set the population size
N_obs_pop = 1000;

% 4.2. Generate data for the independent variable
X_pop = [random('Uniform',-1,1,[N_obs_pop,1])];

% 4.3. Heteroskedasticity parameter 
Gamma = 1.5; 

% 4.4. Generate error
u_pop = random('Normal',0,exp(X_pop*Gamma),[N_obs_pop,1]);

% 4.5. Set true beta of the model
B_true = 0.5; 

% 4.6. Generate y
y_pop = X_pop*B_true+u_pop;

% 4.7. Create the (paired) population data
data_pop = [y_pop X_pop];

%% 5. Plot the scatter diagram and the OLS fitted line
figure;
scatter(X_pop,y_pop,'filled','MarkerFaceColor','black');
hold on
h = lsline;
set(h,'Color','blue','LineWidth',2);
xlabel('X');
ylabel('y');
title('Fig. 1. Simulated heteroskedastic data');
hold off;

%% 6. Draw samples from the population

% 6.1. Set the sample size
N_obs_sample = 100;

% 6.2. Preallocate matrix to save (paired) samples
data_samples_pop = NaN(N_obs_sample,2,N_sim);

% 6.3. Preallocate vector to store coeficient estimates
B_hats_data_samples_pop = NaN(N_sim,1);

% 6.4. Monte Carlo sampling
for i = 1:N_sim
    sample_i = datasample(data_pop,N_obs_sample,'Replace',false);
    data_samples_pop(:,:,i) = sample_i; % Save samples in third dimension
    y = sample_i(:,1);
    X = sample_i(:,2);
    LSS = exercisefunctionlss(y, X);
    B_hats_data_samples_pop(i) = LSS.B_hat(1,1);
end

%% 7. Pick an "initial" sample

% Pick a sample from the samples drawn from the population
data_sample = datasample(data_samples_pop(:,:,i),N_obs_sample, ...
    'Replace',false);

%% 8. Draw (bootsrap) samples from the initial sample 

% 8.1. Preallocate vector to store coefficient estimates
B_hats_data_samples_boot = NaN(N_sim,1);

% 8.2. Draw samples from the original sample and compute the coefficient estimate each time 
for i = 1:N_sim
    data_samples_boot = datasample(data_sample,N_obs_sample,'Replace',true);
    y = data_samples_boot(:,1);
    X = data_samples_boot(:,2);
    LSS = exercisefunctionlss(y,X); 
    B_hats_data_samples_boot(i) = LSS.B_hat(1,1); 
end

%% 9. Plot the estimated PDFs
figure;
hold on
ksdensity(B_hats_data_samples_boot,'function','pdf');
ksdensity(B_hats_data_samples_pop,'function','pdf');
xline(B_true,'Color','black');
ylabel('Density');
xlabel('B\_hat');
legend('Bootstrap sampling','Population sampling','B\_true');
title('Fig. 2. PDF comparison: bootstrap B\_hat vs. population sample B\_hat');
hold off;

%% 10. Plot the estimated CDFs
figure;
hold on
ksdensity(B_hats_data_samples_boot,'function','cdf');
ksdensity(B_hats_data_samples_pop,'function','cdf');
ylabel('Cumulative distribution');
xlabel('B\_hat');
legend('Bootstrap sampling','Population sampling');
title('Fig. 3. CDF comparison: bootstrap B\_hat vs. population sample B\_hat');
hold off
