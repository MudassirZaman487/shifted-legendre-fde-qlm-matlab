function [X,Y] = qlm_second_order_coupled_system(N,sigma1,sigma2,eta, ...
    X0,X1,Y0,Y1,q1,q2,ITT,x)

X = cell(ITT,1);
Y = cell(ITT,1);

a11 = sym(0);
a12 = sym(0);
a21 = sym(0);
a22 = sym(0);

[X{1},Y{1}] = solve_second_order_coupled_linear_legendre( ...
    X0,X1,Y0,Y1,q1,q2,a11,a12,a21,a22,sigma1,sigma2,eta,N,x);

X{1} = expand(vpa(X{1},8));
Y{1} = expand(vpa(Y{1},8));

for i = 1:ITT-1

    Xi = expand(vpa(X{i},8));
    Yi = expand(vpa(Y{i},8));

    A11 = Yi;
    A12 = Xi;

    A21 = Yi;
    A22 = Xi;

    R1 = expand(q1 - Xi*Yi);
    R2 = expand(q2 - Xi*Yi);

    [X{i+1},Y{i+1}] = solve_second_order_coupled_linear_legendre( ...
        X0,X1,Y0,Y1,R1,R2,A11,A12,A21,A22,sigma1,sigma2,eta,N,x);

    X{i+1} = expand(vpa(X{i+1},8));
    Y{i+1} = expand(vpa(Y{i+1},8));

    disp(['QLM iteration = ',num2str(i+1),', N = ',num2str(N), ...
          ', sigma1 = ',num2str(sigma1),', sigma2 = ',num2str(sigma2)])

end

end
