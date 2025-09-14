% Exercise - Understanding the LLN using coin toss simulation

%% 1. Aim of the exercise
% Understanding the law of large numbers using coin toss simulation. 

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Convergence of the sample mean of coin tosses to the population mean

% 3.1. Clear the memory
clear;

% 3.2. Set the number of tosses
N_tosses = 25000; % What would happen if you increase this number? 

% 3.3. Create a sequence of tosses
sequence_of_tosses = random('Binomial',1,0.5,[N_tosses,1]);

% 3.4. Preallocate a matrix
mean_sample = NaN(N_tosses,1); 

% 3.5. Create the sequence of sample means
for i = 1:N_tosses
    mean_sample(i,1) = mean(sequence_of_tosses(1:i,1));
end

% 3.6. Define the population mean
mean_population = 0.5; 

% 3.7. Plot the sample mean and the population mean against the number 
% of tosses
figure
hold on
scatter(1:N_tosses,mean_sample,0.001,'Marker','.')
ylim([0 1]) % Set the range of y axis. 
line([0,N_tosses],[mean_population,mean_population],'Color','red') 
title('Fig. 1. Convergence of the sample mean to the population mean')
legend('Sample means of cumulating tosses',['Population mean ' ...
    'as the sample mean of cumulating tosses in the limit'])
% Or, probability of heads if 1 represents heads. 
ylabel('Sample means of cumulating tosses') 
xlabel('Number of tosses')
hold off

%% 4. Indicator of whether the absolute difference is small enough

% 4.1. Define the absolute difference and the tolerance parameter
abs_dif = abs(mean_sample-mean_population);
e = 0.01;

% 4.2. Create the indicator
abs_dif_smaller_than_e_dummy = abs_dif < e;

% 4.3. Plot the indicator against the number of tosses
figure
hold on
scatter(1:N_tosses,abs_dif_smaller_than_e_dummy,0.001,'Marker','.')
ylim([-1 2]) 
title(['Fig. 2. Cases where the sample mean is close to ' ...
    'the population mean'])
legend('1 if abs\_dif < e, 0 if abs\_dif > e')
ylabel('Whether the absolute difference is smaller than e') 
xlabel('Number of tosses')
hold off

%% 5. Convergence of the sample mean to the population mean in probability

% 5.1. Preallocate the matrix for storing probabilities
prob_of_abs_dif_smaller_than_e = NaN(N_tosses,1);

% 5.2. Probability that the absolute difference is small enough
for i = 1:N_tosses
    prob_of_abs_dif_smaller_than_e(i,1) = mean( ...
        abs_dif_smaller_than_e_dummy(1:i,1)); 
end

% 5.3. Plot the probability that the absolute difference is small enough
figure 
hold on
scatter(1:N_tosses,prob_of_abs_dif_smaller_than_e,0.001,'Marker','.')
ylim([0 1]) 
title(['Fig. 3. Convergence of the sample mean to the pop. mean ' ...
    'in probability'])
ylabel('Probability of abs. dif. smaller than e') 
xlabel('Number of tosses')
hold off
