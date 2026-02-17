clc; close all; clear all;

% -----------------------------
% Basic Parameters
% -----------------------------
fs = 100;
dt = 1/fs;
c = 0.5;
Vsyringe = 3;

T = [2 3 4 5 6];     
n = 5;

% -----------------------------
% Generate Voltage Data
% -----------------------------
for j = 1:n
    t = 0:dt:T(j);
    for i = 1:length(t)
        V(j,i) = c + (1.6-0.2*j)*exp(-t(i));
    end
    L(j) = length(t);
end

% -----------------------------
% Bisection Initialization
% -----------------------------
ml = 0.5; 
mu = 3;
err = 0.001;
iter = 1;

fprintf('\nIter\tb_mid\tf(b)\tError\n');

while abs(mu-ml) > err
    
    mr = (ml+mu)/2;
    
    % ---- Calculate integrals ----
    for j = 1:n
        I(j) = 0;
        for i = 1:L(j)
            I(j) = I(j) + (V(j,i)-c)^mr * dt;
        end
    end
    
    % ---- Standard Deviation ----
    meanI = sum(I)/n;
    f = sqrt(sum((I-meanI).^2)/n);
    
    % ---- Store error for graph ----
    error_list(iter) = abs(mu-ml);
    b_list(iter) = mr;
    
    % ---- Bisection update ----
    if f > 0
        mu = mr;
    else
        ml = mr;
    end
    
    fprintf('%d\t%.4f\t%.6f\t%.6f\n',iter,mr,f,abs(mu-ml));
    
    iter = iter + 1;
end

b_final = mr;

% -----------------------------
% Scaling Factor A
% -----------------------------
A = Vsyringe / I(1);

fprintf('\nFinal b = %.6f\n',b_final);
fprintf('Scaling A = %.6f\n',A);

% -----------------------------
% Plot 1: Objective Function
% -----------------------------
b_vals = 0.5:0.05:3;

for k = 1:length(b_vals)
    for j = 1:n
        Itemp(j) = 0;
        for i = 1:L(j)
            Itemp(j) = Itemp(j) + (V(j,i)-c)^b_vals(k)*dt;
        end
    end
    mtemp = sum(Itemp)/n;
    f_vals(k) = sqrt(sum((Itemp-mtemp).^2)/n);
end

figure;
plot(b_vals,f_vals,'b','LineWidth',2);
hold on;
plot(b_final,0,'ro','MarkerSize',8,'LineWidth',2);
xlabel('Exponent b');
ylabel('f(b)');
title('Objective Function vs Exponent b');
legend('f(b)','Intersection (Optimal b)','Location','best');
grid on;

% -----------------------------
% Plot 2: Error Convergence
% -----------------------------
figure;
plot(1:length(error_list),error_list,'r','LineWidth',2);
hold on;

% Final iteration point (intersection with tolerance)
final_iter = length(error_list);
final_error = error_list(end);

plot(final_iter,final_error,'bo','MarkerSize',8,'LineWidth',2);

% Tolerance line
yline(err,'k--','Tolerance');

xlabel('Iteration');
ylabel('Error |b_u - b_l|');
title('Error Convergence of Bisection Method');
legend('Error','Final Converged Point','Tolerance','Location','best');
grid on;