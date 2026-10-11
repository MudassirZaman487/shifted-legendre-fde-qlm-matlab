clc
clear all
close all
syms x
eta = 1;
f = cos(2*pi*x);
for N = 4:2:12
    A = shifted_legendre_fun_vector(N,eta,x);
    C = shifted_legendre_coef(f,N,eta,x);
    f_app = C*A;
figure(1)
subplot(1,2,1)
ezplot(f_app,[0,eta])
hold on
grid on
title('Exact and approximate solutions')
xlabel('x')
ylabel('f(x)')
axis auto
subplot(1,2,2)
hold on
grid on
ezplot(abs(f-f_app),[0,eta])
set(gca,'YScale','Log')
title('Absolute error')
xlabel('x')
ylabel('|f(x)-f_N(x)|')
axis auto
N
pause(1)
end

dom = linspace(0,eta,10);
subplot(1,2,1)
hold on
plot(dom,subs(f,x,dom))
