clc
clear all
close all

syms x
eta = 1;
zeta = 1;

for N = 5:10
    u = sin(x);
    g = x^3*(1-x);
    Exact = g*diff(u,x);
    A = shifted_legendre_fun_vector(N,eta,x);
    U = double(shifted_legendre_coef(u,N,eta,x));
    DiffMat = double(shifted_legendre_diff(N,zeta,eta));
    VarMat = double(shifted_legendre_var_matrix(g,N,eta,x));
    App = U*DiffMat*VarMat*A;
    figure(1)
    hold on
    ezplot(abs(Exact-App),[0,eta])
    grid on
    title('(a)')
    xlabel('x')
    ylabel('Absolute error')
    axis auto
    box on
    N
    pause(1)
end
