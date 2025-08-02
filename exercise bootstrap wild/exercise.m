% Exercise - Understanding the method of wild bootstrap

%% 1. Aim of the Exercise
% To learn about the method of wild bootstrap.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Set values for the parameters of the simulation

% 3.1. Clear the memory
clear;

% 3.2. Set the number of simulations as the number of bootstrap samples
N_sim = 5000;

%% 4. Generate population data

% 4.1. Set the population size
N_obs_pop = 10000;

% 4.2. Generate data for the independent variable
X_pop = random('Uniform',-1,1,[N_obs_pop 1]);

% 4.3. Heteroskedasticity parameter 
Gamma = 1.5; 

% 4.4. Generate error
u_pop = random('Normal',0,exp(X_pop*Gamma),[N_obs_pop 1]);

% 4.5. Set true beta of the model
B_true = 0.5; 

% 4.6. Generate y
y_pop = X_pop*B_true+u_pop;

% 4.7. Generate (paired) population data
data_pop = [y_pop X_pop];

%% 5. Plot the scatter diagram and the OLS fitted line
figure;
scatter(X_pop,y_pop,'filled','MarkerFaceColor','black');
hold on
set(lsline,'color','blue','LineWidth',2);
xlabel('X');
ylabel('y');
title('Fig. 1. Simulated heteroskedastic data');
hold off;

%% 6. Generate a sample

% 6.1. Set the sample size
N_obs_sample = 500;

% 6.2. Draw a sample from the population
data_sample = datasample(data_pop,N_obs_sample,'Replace',false);

%% 7. Draw (bootstrap) samples from the sample 

% 7.1. Create y and X
y = data_sample(:,1);
X = data_sample(:,2);

% 7.2. Estimate beta
LSS = exercisefunctionlss(y,X); 
B_hat = LSS.B_hat(1,1);

% 7.3. Compute the residuals
u_hat = LSS.u_hat;

% 7.4. Compute leverage values used in HC2 correction
h = diag(X/(X'*X)*X');

% 7.5. Preallocate vector to store coefficient estimates
B_hats_data_samples_boot = NaN(N_sim,1);

% 7.6. Draw samples using residual multipliers and compute the coefficient estimate each time
for i = 1:N_sim
    v = random('Normal',0,1,[N_obs_sample 1]); % N(0,1) multipliers
    % v = 2*(rand(N_obs_sample,1) > 0.5)-1; % Rademacher multipliers
    y_boot = X*B_hat+(u_hat./sqrt(1-h)).*v; % u_hat./sqrt(1-h) is the HC2 adjustment applied to residuals
    LSS = exercisefunctionlss(y_boot,X); 
    B_hats_data_samples_boot(i) = LSS.B_hat(1,1);
end

%% 8. Draw samples from the population

% 8.1. Preallocate vector to store coeficient estimatess
B_hats_data_samples_pop = NaN(N_sim,1);

% 8.2. Draw samples from the population and compute the coefficient estimate each time
for i = 1:N_sim
    data_sample = datasample(data_pop,N_obs_sample,'Replace',false);
    y = data_sample(:,1);
    X = data_sample(:,2);
    LSS = exercisefunctionlss(y,X); 
    B_hats_data_samples_pop(i) = LSS.B_hat(1,1);
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
hold off

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
