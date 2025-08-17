% Exercise - Halton sequences

%% 1. Aim of the exercise 
% Halton sequences are low-discrepancy, quasi-random sequences used for
% uniform sampling in high-dimensional spaces. This exercise illustrates
% their structure and behavior through visual comparisons with uniform and
% normal distributions, using both 2D and 3D plots.

%% 2. Theory
% Refer to the accompanying PDF file for the theory.

%% 3. Set the paramaters of the Halton sequence

% 3.1. Clear the memory 
clear;

% 3.2. Set the number of observational units
N = 100;

% 3.3. Set the number of dimensions
dimensions = 1;

% 3.4. Set the number of draws per observational unit
draws = 1;

% 3.5. Define prime base for Halton sequence generation
prnum = 3;

%% 4. Sampling with Halton sequences and uniform distribution

% 4.1. Generate a Halton sequence
[H,~] = halton(N,dimensions,draws,'prnum',prnum);

% 4.2. Generate samples from a uniform(0,1) distribution
U = random("Uniform",0,1,[N,1]);

%% 5. Plot the Halton and uniform draws

% 5.1. Plot
figure
hold on
plot(H,0.4*ones(N,1),'ro','MarkerFaceColor','red');
plot(U,0.6*ones(N,1),'go','MarkerFaceColor','green');
axis([0 1 0 1])
legend("Halton draws","Uniform draws")
title("Fig. 1. Comparison of Halton and uniform draws")
hold off

%% 6. Set the paramaters of the Halton sequences

% 6.1. Clear the memory 
clear;

% 6.2. Set the number of observational units
N = 100;

% 6.3. Set the number of dimensions
dimensions = 2;

% 6.4. Set the number of draws per observational unit
draws = 1;

% 6.5. Define prime bases for Halton sequence generation
prnum = [43,47];

%% 7. Halton sequences and transformation to standard normal distribution

% 7.1. Generate regular Halton sequence and transform to standard normal
[H_standard_2D,Z_standard_2D] = halton(N,dimensions,draws, ...
    'prnum',prnum);

% 7.2. Generate scrambled Halton sequence and transform to standard normal
[H_scramble_2D,Z_scramble_2D] = halton(N,dimensions,draws, ...
    'prnum',prnum,'scramble',1);

% 7.3. Generate leap Halton sequence (every 5th entry) and transform to standard normal
[H_leap_2D,Z_leap_2D] = halton(N,dimensions,draws, ...
    'prnum',prnum,'leap',4);

% 7.4. Generate randomized Halton sequence and transform to standard normal
[H_random_2D,Z_random_2D] = halton(N,dimensions,draws, ...
    'prnum',prnum,'random',1);

%% 8. Visual comparison of standard and scrambled Halton draws

% 8.1. Plot standard Halton sequence in 2D
figure
plot(H_standard_2D(:,1),H_standard_2D(:,2),'ro','MarkerFaceColor','red')
xlabel("Dimension 1")
ylabel("Dimension 2")
legend("Standard Halton sequence")
title("Fig. 2. Standard Halton draws")

% 8.2. Plot scrambled Halton sequence in 2D
figure
plot(H_scramble_2D(:,1),H_scramble_2D(:,2),'go','MarkerFaceColor','green')
xlabel("Dimension 1")
ylabel("Dimension 2")
legend("Scrambled Halton sequence")
title("Fig. 3. Scrambled Halton draws")

%% 9. Density estimation and 3D surface plot setup

% 9.1. Define grid for evaluating density surfaces
x1 = linspace(min(Z_standard_2D(:,1)),max(Z_standard_2D(:,1)),50);
x2 = linspace(min(Z_standard_2D(:,2)),max(Z_standard_2D(:,2)),50);
[X1,X2] = meshgrid(x1,x2);
grid_points = [X1(:),X2(:)];

% 9.2. Turn the transformed Halton draws into a density function
[f_standard,~] = ksdensity(Z_standard_2D,grid_points,'Function','pdf');
[f_scramble,~] = ksdensity(Z_scramble_2D,grid_points,'Function','pdf');
[f_leap,~] = ksdensity(Z_leap_2D,grid_points,'Function','pdf');
[f_random,~] = ksdensity(Z_random_2D,grid_points,'Function','pdf');

% 9.3. Reshape the densities for the surface plot
F_standard = reshape(f_standard,size(X1));
F_scrambled = reshape(f_scramble,size(X1));
F_leap = reshape(f_leap,size(X1));
F_random = reshape(f_random,size(X1));

% 9.4. Compute reference density from bivariate standard normal distribution
Z_reference = mvnpdf(grid_points);
Z_reference = reshape(Z_reference,size(X1));

%% 10. Compare Halton-based sampling densities to multivariate normal in 3D

% 10.1. Standard Halton vs. multivariate normal density
figure
hold on
surf(X1,X2,F_standard,'FaceColor','y','FaceAlpha',0.7);
surf(X1,X2,Z_reference,'FaceColor','m','FaceAlpha',0.1);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Density');
legend("Standard Halton density","Multivariate normal PDF");
title("Fig. 4. Standard Halton sampling vs. normal density");
view(3);
grid on
hold off

% 10.2. Scrambled Halton vs. multivariate normal density
figure
hold on
surf(X1,X2,F_scrambled,'FaceColor','y','FaceAlpha',0.7);
surf(X1,X2,Z_reference,'FaceColor','m','FaceAlpha',0.1);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Density');
legend("Scrambled Halton density","Multivariate normal PDF");
title("Fig. 5. Scrambled Halton sampling vs. normal density");
view(3);
grid on
hold off

% 10.3. Leap Halton vs. multivariate normal density
figure
hold on
surf(X1,X2,F_leap,'FaceColor','y','FaceAlpha',0.7);
surf(X1,X2,Z_reference,'FaceColor','m','FaceAlpha',0.1);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Density');
legend("Leap Halton density","Multivariate normal PDF");
title("Fig. 6. Leap Halton sampling vs. normal density");
view(3);
grid on
hold off

% 10.4. Randomized Halton vs. multivariate normal density
figure
hold on
surf(X1,X2,F_random,'FaceColor','y','FaceAlpha',0.7);
surf(X1,X2,Z_reference,'FaceColor','m','FaceAlpha',0.1);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Density');
legend("Randomized Halton density","Multivariate normal PDF");
title("Fig. 7. Randomized Halton sampling vs. normal density");
view(3);
grid on
hold off
