clc
clear all
close all
format long

syms x

eta = 1;
N_fixed = 10;
ITT = 8;

N_values = 5:1:12;
gamma_values = 0.5:0.1:1;
dom = linspace(0,eta,400);

if ~exist('Thesis Graphs','dir')
    mkdir('Thesis Graphs');
end
u_exact = (exp(2*x)-1)/(exp(2*x)+1);
u0 = 0;
u_exact_val = double(subs(u_exact,x,dom));

%% ==========================================================
%  Figure 1
%  Approaching behavior of QLM iterations to exact solution
%  Fixed scale level N = 10 and gamma = 1
% ==========================================================
gamma = 1;
U = qlm_riccati_first_order(N_fixed,gamma,eta,u0,ITT,x);
rmse_iteration = zeros(ITT,1);
figure(1)
subplot(1,3,1)
hold on
grid on
box on
for i = 1:ITT
    Ui_val = double(subs(U{i},x,dom));
    ezplot(U{i},[0,eta])
end
dom_ex = linspace(0,eta,10);
ezplot(u_exact,[0,eta])
title('(a)')
xlabel('x')
ylabel('u(x)')
legend_text = cell(1,ITT+1);
for i = 1:ITT
    legend_text{i} = ['Iter = ',num2str(i)];
end
legend_text{ITT+1} = 'Exact';
legend(legend_text,'Location','Best')
subplot(1,3,2)
hold on
grid on
box on

for i = 1:ITT
    Ui_val = double(subs(U{i},x,dom));
    rel_err = abs(u_exact-U{i})/(abs(u_exact));
    ezplot(rel_err,[0,eta])
    rmse_iteration(i) = sqrt(mean((u_exact_val-Ui_val).^2));
end

title('(b)')
xlabel('x')
ylabel('Relative absolute error')

legend_text = cell(1,ITT);
for i = 1:ITT
    legend_text{i} = ['Iter = ',num2str(i)];
end
legend(legend_text,'Location','Best')

subplot(1,3,3)
semilogy(1:ITT,rmse_iteration,'o-','LineWidth',1.5)
grid on
box on
title('(c)')
xlabel('Iteration number')
ylabel('RMSE')

%% ==========================================================
%  Figure 2
%  Convergence for different scale levels N
%  Fixed iteration number ITT
% ==========================================================

rmse_N = zeros(length(N_values),1);

ITT = 4;
figure(2)

subplot(1,2,1)
hold on
grid on
box on

for i = 1:length(N_values)
    N = N_values(i);
    U_N = qlm_riccati_first_order(N,gamma,eta,u0,ITT,x);
    U_final = U_N{ITT};
    U_val = double(subs(U_final,x,dom));
    rel_err = abs(u_exact-U_final)/(abs(u_exact));
    rmse_N(i) = sqrt(mean((u_exact_val-U_val).^2));
    ezplot(rel_err,[0,eta])
    set(gca,'YScale', 'Log')
    axis auto
i
pause(1)
end
title('(a)')
xlabel('x')
set(gca,'YScale', 'Log')
ylabel('Relative absolute error')

legend_text = cell(1,length(N_values));
for i = 1:length(N_values)
    legend_text{i} = ['N = ',num2str(N_values(i))];
end
legend(legend_text,'Location','Best')

subplot(1,2,2)
semilogy(N_values,rmse_N,'o-','LineWidth',1.5)
grid on
box on
title('(b)')
xlabel('Scale level N')
ylabel('RMSE')


% % ==========================================================
%  Figure 3
%  Fractional-order approaching behavior as gamma approaches 1
% ==========================================================
% 
figure(3)
hold on
grid on
box on
N = 9;
ITT = 4;
for i = 1:length(gamma_values)
    gamma = gamma_values(i);
    U_gamma = qlm_riccati_first_order(N,gamma,eta,u0,ITT,x);
    U_final = U_gamma{ITT};
    ezplot(U_final,[0,eta])
i
pause(1)
end
dom_ex = linspace(0,eta,10);
plot(dom_ex, subs(u_exact,x,dom_ex))
title('Fractional-order approaching behavior')
xlabel('x')
ylabel('u(x)')
legend_text = cell(1,length(gamma_values)+1);
for i = 1:length(gamma_values)
    legend_text{i} = ['\gamma = ',num2str(gamma_values(i),'%.1f')];
end
legend_text{length(gamma_values)+1} = 'Exact at \gamma = 1';
legend(legend_text,'Location','Best')

