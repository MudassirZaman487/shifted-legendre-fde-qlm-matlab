clc
clear all
close all

syms x

N = 8;
eta = 1;
zeta_val = 0.5:0.1:1;

f = x*(x-0.5)*(x-1);
Fexact = int(f,x,0,x)

for i = 1:length(zeta_val)
    zeta = zeta_val(i);
    A = shifted_legendre_fun_vector(N,eta,x);
    C = double(shifted_legendre_coef(f,N,eta,x));
    IntMat = double(shifted_legendre_int(N,zeta,eta));
    Fmat = C*IntMat*A;
    figure(1)
    hold on
    ezplot(Fmat,[0,eta])
    grid on
    title('(a)')
    xlabel('x')
    ylabel('Fractional integral')
    axis auto
    box on
    zeta
    pause(1)
end

dom = linspace(0,eta,10);

figure(1)
hold on
plot(dom,subs(Fexact,x,dom))
