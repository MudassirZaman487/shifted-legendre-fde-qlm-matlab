clc
clear all
close all
format long

syms x

eta = 1;
N_fixed = 8;
ITT = 8;

rho_values = 0.5:0.1:1;


dom = linspace(0,eta,400);

if ~exist('Thesis Graphs','dir')
    mkdir('Thesis Graphs');
end

U_exact = 8*x^3*(1-x)+(1-x)^2+3;

a_x = x;
b_x = x^2;

Ux = diff(U_exact,x,1);
Uxx = diff(U_exact,x,2);

f = Uxx - (b_x*Ux + a_x*U_exact + U_exact^2);

y0 = subs(U_exact,x,0);
y1 = subs(diff(U_exact,x,1),x,0);

U_exact_val = double(subs(U_exact,x,dom));

%% ==========================================================
%  Figure 1
%  Approaching behavior of QLM iterations to exact solution
%  Fixed scale level N = 8, sigma = 2, rho = 1
% ==========================================================

sigma = 2;
rho = 1;

U = qlm_second_order_single(N_fixed,sigma,rho,eta,y0,y1,f,a_x,b_x,ITT,x);

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
dom_ex= linspace(0,eta,10);
plot(dom_ex, subs(U_exact,x,dom_ex))

title('(a)')
xlabel('x')
ylabel('U(x)')

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
    rel_err = abs(U_exact-U{i})/(abs(U_exact));
    ezplot(rel_err,[0,eta])
    rmse_iteration(i) = sqrt(mean((U_exact_val-Ui_val).^2));
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

subplot(1,3,3)
semilogy(1:ITT,rmse_iteration,'o-','LineWidth',1.5)
grid on
box on
title('(c)')
xlabel('Iteration number')
ylabel('RMSE')



%% ==========================================================
%  Figure 2
%  Approaching behavior as rho approaches 1
%  Fixed sigma = 2
% ==========================================================

figure(3)
hold on
grid on
box on

N = 5;
ITT = 5;
sigma = 2;

for i = 1:length(rho_values)

    rho = rho_values(i);

    U_rho = qlm_second_order_single(N,sigma,rho,eta,y0,y1,f,a_x,b_x,ITT,x);
    U_final = U_rho{ITT};

    ezplot(U_final,[0,eta])

    i
    pause(1)

end

dom_ex= linspace(0,eta,10);
plot(dom_ex, subs(U_exact,x,dom_ex));
title('Approaching behavior as \rho approaches 1')
xlabel('x')
ylabel('U(x)')

legend_text = cell(1,length(rho_values)+1);
for i = 1:length(rho_values)
    legend_text{i} = ['\rho = ',num2str(rho_values(i),'%.1f')];
end
legend_text{length(rho_values)+1} = 'Exact at \rho = 1';
legend(legend_text,'Location','Best')




