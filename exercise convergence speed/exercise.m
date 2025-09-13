% Exercise - Understanding convergence speed through Big O notation

%% 1. Aim of the exercise
% This exercise demonstrates how convergence rates depend on dimension,
% using the expression 1/n^(1/d) as a simplified model of error reduction.
% As the number of samples n increases, error decreases, but this happens
% more slowly in higher dimensions. This effect is known as the curse of
% dimensionality. We compare convergence behavior for dimensions d = 1, 2,
% and 3, and visualize how the rate of decay changes. The goal is to build
% intuition for how Big O notation reflects convergence speed in numerical
% analysis.

%% 2. Define sample size range

% Sample size range
n = 1:1:100; % Avoid n = 0 to prevent division by zero.

%% 3. Model convergence rates

% Convergence rates for d = 1, 2, 3
rate_d1 = 1./(n.^(1/1)); 
rate_d2 = 1./(n.^(1/2));
rate_d3 = 1./(n.^(1/3));

%% 4. Plot convergence behavior across dimensions

% Plot
figure
set(gcf,'Position',[100,100,1000,1000]); 

hold on
plot(n,rate_d1,'DisplayName','d = 1');
plot(n,rate_d2,'DisplayName','d = 2');
plot(n,rate_d3,'DisplayName','d = 3');
title('Fig. 1. Convergence rates modeled by 1/n^{1/d}');
xlabel('Sample size (n)');
ylabel('Estimated error');
legend('show');
hold off
