% Exercise - Understanding the theory of importance sampling

%% 1. Aim of the exercise

% The aim of this exercise is to understand the role of tail behaviour in
% importance sampling, by comparing proposal and target distributions and
% observing how heavier tails help ensure finite variance of the estimator.

%% 2. Tails comparison

% 2.1. Sequence from -5 to 5 by steps of 0.1
x_values = -5:0.1:5; 

% 2.2. Theoretical Standard Normal PDF
PD_Normal = makedist('Normal','mu',0,'sigma',1);

% 2.3. Calculate the PDF for the Normal distribution
PDF_Normal = pdf(PD_Normal,x_values);

% 2.4. Theoretical standard cauchy PDF
PD_Cauchy = makedist('tLocationScale','mu',0,'sigma',1,'nu',1);

% 2.5. Title
PDF_Cauchy = pdf(PD_Cauchy,x_values);

%% 3. Plot the theoretical PDFs

% Plot
figure
hold on
plot(x_values,PDF_Normal,'DisplayName','N(0,1)')
plot(x_values,PDF_Cauchy,'DisplayName','Cauchy(0,1)')
title('Fig. 1. Importance sampling PDF tails')
xlabel('x')
ylabel('PDF')
legend('show');
hold off
