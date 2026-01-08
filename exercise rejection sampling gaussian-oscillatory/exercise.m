% Exercise - Understanding rejection sampling using a Gaussian-oscillatory

%% 1. Aim of the exercise

% The aim of this exercise is to develop intuition for rejection sampling.
% Specifically, we use a Gaussian envelope combined with an oscillatory
% component as the target function, and apply rejection sampling to
% estimate the area under this function. By comparing accepted and rejected
% samples under different proposal distributions, the exercise illustrates
% how the choice of proposal affects efficiency and accuracy.

%% 2. Theory

% Refer to the accompanying PDF file for the theory.

%% 3. Clear workspace variables

% Clear the memory 
clear;

%% 4. Define a domain of values

% 4.1. Specify number of random samples to estimate the integral
N_samples = 1000;

% 4.2. Define a domain of values
x = -4:0.001:4;

%% 5. Construct the target function

% 5.1. Define the oscillatory component  
oscillatory_component = pi*cos(x).^2.*sin(pi*x).^2 + 1;

% 5.2. Preallocate the target function vector
target_function = NaN(1,length(x));

% 5.3. Define the target function piecewise  
for i = 1:length(x)
    % Check if x(i) lies within the domain [-3, 3]
    if abs(x(i)) <= 3
        % Inside domain: compute Gaussian * Oscillatory component
        target_function(i) = exp(-0.5*x(i).^2).* ...
            oscillatory_component(i);
    else
        % Outside domain: set target function to zero
        target_function(i) = 0;
    end
end

%% 6. Construct a proposal distribution

% 6.1. Specify the proposal distribution  
proposal_uniform_PD = makedist('Uniform',-3,3);

% 6.2. Compute the proposal PDF  
proposal_uniform_PDF = pdf(proposal_uniform_PD,x);

% 6.3. Determine the maximum of the target function  
c = max(target_function);

% 6.4. Choose envelope constant 
C = c*(3-(-3));

%% 7. Visualize the target and proposal distributions

% Create new figure
figure
hold on
% Plot the target function (blue curve)
plot(x,target_function, ...
    'Color','blue', ...
    'DisplayName','Target function')
plot(x,proposal_uniform_PDF, ...
    'DisplayName','Proposal PDF')
% Plot the scaled proposal distribution (black curve)
plot(x,C*proposal_uniform_PDF, ...
    'Color','black', ...
    'DisplayName','Scaled proposal PDF')
plot(-4:-3,[max(target_function) max(target_function)],'--', ...
    'DisplayName','Max target function')
xline(-1.2,'--', ...
    'DisplayName','Reference x = -1.2')
f_Value = target_function((-1.2+4)/0.01); % Intersection point f(-1.2)
plot(-1.2,f_Value,'Marker','.', ...
    'DisplayName','Intersection f(-1.2)')
text(-5.1,max(target_function),'C \cdot p_X(-1.2)')
text(-1.3,-0.08,'-1.2')
title('Fig. 1. Target vs Proposal Distributions in Rejection Sampling')
xlabel('x')
legend('show');
hold off

%% 8. Determine accepted and rejected samples

% 8.1. Generate x-coordinates from the uniform proposal on [-3,3]
proposal_samples = random('Uniform',-3,3,[N_samples 1]);

% 8.2. Generate U ~ Unif(0,1) for vertical placement
u_samples = random('Uniform',0,1,[N_samples 1]);

% 8.3. Compute corresponding y-coordinates scaled by envelope height
y_coordinates = u_samples.* ...
    (C*pdf(proposal_uniform_PD,proposal_samples));

% 8.4. Combine x- and y-coordinates
proposal_points_scaled = [proposal_samples y_coordinates];

% 8.5. Preallocate matrices for accepted and rejected points
accepted_samples = NaN(N_samples,2);
rejected_samples = NaN(N_samples,2);

% 8.6. Classify each sample
for j = 1:N_samples
    x_coordinate = proposal_points_scaled(j,1);
    y_coordinate = proposal_points_scaled(j,2);
    % Find index in x-values for target_function(x_coordinate)
    find_index = (abs(min(x)) + x_coordinate) / 0.001;
    index = int16(find_index) + 1; % Adjust for round-off error
    if y_coordinate <= target_function(index)
        % Accepted sample
        accepted_samples(j,:) = [x_coordinate, y_coordinate];
    else
        % Rejected sample
        rejected_samples(j,:) = [x_coordinate, y_coordinate];
    end
end

%% 9. Visualize accepted and rejected samples

% Create new figure
figure
hold on
% Plot the target function (blue curve)
plot(x,target_function, ...
    'Color','blue', ...
    'DisplayName','Target function')
% Plot the scaled proposal distribution (black curve)
plot(x,C*proposal_uniform_PDF, ...
    'Color','black', ...
    'DisplayName','Scaled proposal')
% Plot accepted points (below target function, blue dots)
scatter(accepted_samples(:,1),accepted_samples(:,2),10, ...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue', ...
    'DisplayName','Accepted samples')
% Plot rejected points (above target function, black dots)
scatter(rejected_samples(:,1),rejected_samples(:,2),10, ...
    'MarkerEdgeColor','black','MarkerFaceColor','black', ...
    'DisplayName','Rejected samples')
title('Fig. 2. Accepted vs Rejected Samples in Rejection Sampling')
xlabel('x')
legend('show')
hold off

%% 10. Construct a new proposal distribution

% 10.1. Define the proposal distribution as a standard normal N(0,1)
proposal_standard_normal_PD = makedist('Normal','mu',0,'sigma',1);

% 10.2. Compute the proposal PDF values over the domain x
proposal_standard_normal_PDF = pdf(proposal_standard_normal_PD,x);

% 10.3. Choose envelope constant C_N so that C_N * proposal_PDF >= 
% target_function
C_N = 9;

%% 11. Visualize the target and new proposal distributions

% Create new figure
figure
hold on
% Plot the target function (blue curve)
plot(x,target_function, ...
    'Color','blue', ...
    'DisplayName','Target function')
% Plot the scaled proposal distribution (black curve)
plot(x,C_N*proposal_standard_normal_PDF, ...
    'Color','black', ...
    'DisplayName','Scaled proposal')
title('Fig. 3. Target vs Scaled Normal Proposal in Rejection Sampling')
xlabel('x')
legend('show')
hold off

%% 12. Determine accepted and rejected samples

% 12.1. Generate x-coordinates (realizations from N(0,1))
proposal_samples = random('Normal',0,1,[N_samples 1]);  

% 12.2. Compute corresponding y-coordinates
proposal_pdf_values = (1/sqrt(2*pi)) * ...
    exp(-0.5 * proposal_samples.^2);

% 12.3. Generate uniform scaling factors to spread points vertically
u_samples = random('Uniform',0,1,[N_samples 1]);

% 12.4. Define envelope constant (from Section 10)
c = C_N;

% 12.5. Scale y-coordinates by envelope constant and uniform factor
proposal_points_scaled = [proposal_samples ...
    c * u_samples .* proposal_pdf_values];

% 12.6. Preallocate matrix for accepted points (below target function)
accepted_samples = NaN(N_samples,2);

% 12.7. Preallocate matrix for rejected points (above target function)
rejected_samples = NaN(N_samples,2);

% 12.8. Loop through all samples and classify each point
for k = 1:N_samples
    x_coordinate = proposal_points_scaled(k,1);
    y_coordinate = proposal_points_scaled(k,2);
    find_index = (abs(min(x)) + x_coordinate) / 0.001; 
    index = int16(find_index) + 1; % Adjust for round-off error
    if y_coordinate <= target_function(index)
        accepted_samples(k,:) = [x_coordinate, y_coordinate];
    else
        rejected_samples(k,:) = [x_coordinate, y_coordinate];
    end
end

%% 13. Visualize accepted and rejected samples

% 13.1. Create a new figure window
figure
hold on
% Plot the target function (blue curve)
plot(x,target_function, ...
    'Color','blue', ...
    'DisplayName','Target function')
% Plot the scaled standard normal proposal distribution (black curve)
plot(x,C_N*proposal_standard_normal_PDF, ...
    'Color','black', ...
    'DisplayName','Scaled proposal')
% Plot accepted samples (blue dots below target function)
scatter(accepted_samples(:,1),accepted_samples(:,2),10, ...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue', ...
    'DisplayName','Accepted samples')
% Plot rejected samples (black dots above target function)
scatter(rejected_samples(:,1),rejected_samples(:,2),10, ...
    'MarkerEdgeColor','black','MarkerFaceColor','black', ...
    'DisplayName','Rejected samples')
title('Fig. 4. Accepted vs Rejected Samples with Normal Proposal')
xlabel('x')
legend('show')
hold off
