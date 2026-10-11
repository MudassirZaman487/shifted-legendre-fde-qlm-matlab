clc
clear all
close all
format long

syms x

eta = 1;
N_fixed = 8;
ITT = 7; % For Figure 1

gamma_values = 0.5:0.1:1;
dom = linspace(0,eta,400);

if ~exist('Thesis Graphs','dir')
    mkdir('Thesis Graphs');
end

X_exact = x*(x-0.5)*(x-1);
Y_exact = x*(x-1)*(x+1);

X0 = subs(X_exact,x,0);
Y0 = subs(Y_exact,x,0);

q1 = diff(X_exact,x,1) - X_exact*Y_exact;
q2 = diff(Y_exact,x,1) - X_exact^2;

X_exact_val = double(subs(X_exact,x,dom));
Y_exact_val = double(subs(Y_exact,x,dom));

%% ==========================================================
%  Figure 1
%  Approaching behavior of QLM iterations to exact solution
%  Fixed scale level N = 8 and gamma = 1
% ==========================================================

mu1 = 1;
mu2 = 1;

[X,Y] = qlm_first_order_coupled_system(N_fixed,mu1,mu2,eta,X0,Y0,q1,q2,ITT,x);

rmse_X_iteration = zeros(ITT,1);
rmse_Y_iteration = zeros(ITT,1);

figure(1)

subplot(2,3,1)
hold on
grid on
box on

for i = 1:ITT
    Xi_val = double(subs(X{i},x,dom));
    ezplot(X{i},[0,eta])
end

dom_ex = linspace(0,eta,10);
hold on
plot(dom_ex, subs(X_exact,x,dom_ex))

title('(a)')
xlabel('x')
ylabel('X(x)')

legend_text = cell(1,ITT+1);
for i = 1:ITT
    legend_text{i} = ['Iter = ',num2str(i)];
end
legend_text{ITT+1} = 'Exact';
legend(legend_text,'Location','Best')

subplot(2,3,2)
hold on
grid on
box on

for i = 1:ITT
    Xi_val = double(subs(X{i},x,dom));
    rel_err_X = abs(X_exact-X{i})/(abs(X_exact));
    ezplot(rel_err_X,[0,eta])
    rmse_X_iteration(i) = sqrt(mean((X_exact_val-Xi_val).^2));
end

title('(b)')
xlabel('x')
ylabel('Relative absolute error')
set(gca,'YScale','Log')

legend_text = cell(1,ITT);
for i = 1:ITT
    legend_text{i} = ['Iter = ',num2str(i)];
end
legend(legend_text,'Location','Best')

subplot(2,3,3)
semilogy(1:ITT,rmse_X_iteration,'o-','LineWidth',1.5)
grid on
box on
title('(c)')
xlabel('Iteration number')
ylabel('RMSE')

subplot(2,3,4)
hold on
grid on
box on

for i = 1:ITT
    Yi_val = double(subs(Y{i},x,dom));
    ezplot(Y{i},[0,eta])
end
dom_ex = linspace(0,eta,10);
plot(dom_ex, subs(Y_exact,x,dom_ex))

title('(d)')
xlabel('x')
ylabel('Y(x)')

legend_text = cell(1,ITT+1);
for i = 1:ITT
    legend_text{i} = ['Iter = ',num2str(i)];
end
legend_text{ITT+1} = 'Exact';
legend(legend_text,'Location','Best')

subplot(2,3,5)
hold on
grid on
box on
for i = 1:ITT
    Yi_val = double(subs(Y{i},x,dom));
    rel_err_Y = abs(Y_exact-Y{i})/(abs(Y_exact));
    ezplot(rel_err_Y,[0,eta])
    rmse_Y_iteration(i) = sqrt(mean((Y_exact_val-Yi_val).^2));
end

title('(e)')
xlabel('x')
ylabel('Relative absolute error')
set(gca,'YScale','Log')

legend_text = cell(1,ITT);
for i = 1:ITT
    legend_text{i} = ['Iter = ',num2str(i)];
end
legend(legend_text,'Location','Best')

subplot(2,3,6)
semilogy(1:ITT,rmse_Y_iteration,'o-','LineWidth',1.5)
grid on
box on
title('(f)')
xlabel('Iteration number')
ylabel('RMSE')




%% ==========================================================
%  Figure 2
%  Fractional-order approaching behavior as gamma approaches 1
% ==========================================================

figure(3)

N = 8;
ITT = 5;

subplot(1,2,1)
hold on
grid on
box on

for i = 1:length(gamma_values)

    gamma = gamma_values(i);

    [X_gamma,Y_gamma] = qlm_first_order_coupled_system(N,gamma,gamma,eta,X0,Y0,q1,q2,ITT,x);

    X_final = X_gamma{ITT};

    ezplot(X_final,[0,eta])

    i
    pause(1)

end

ezplot(X_exact,[0,eta])

title('(a)')
xlabel('x')
ylabel('X(x)')

legend_text = cell(1,length(gamma_values)+1);
for i = 1:length(gamma_values)
    legend_text{i} = ['\gamma = ',num2str(gamma_values(i),'%.1f')];
end
legend_text{length(gamma_values)+1} = 'Exact at \gamma = 1';
legend(legend_text,'Location','Best')

subplot(1,2,2)
hold on
grid on
box on

for i = 1:length(gamma_values)

    gamma = gamma_values(i);

    [X_gamma,Y_gamma] = qlm_first_order_coupled_system(N,gamma,gamma,eta,X0,Y0,q1,q2,ITT,x);

    Y_final = Y_gamma{ITT};

    ezplot(Y_final,[0,eta])

    i
    pause(1)

end

ezplot(Y_exact,[0,eta])

title('(b)')
xlabel('x')
ylabel('Y(x)')

legend_text = cell(1,length(gamma_values)+1);
for i = 1:length(gamma_values)
    legend_text{i} = ['\gamma = ',num2str(gamma_values(i),'%.1f')];
end
legend_text{length(gamma_values)+1} = 'Exact at \gamma = 1';
legend(legend_text,'Location','Best')

