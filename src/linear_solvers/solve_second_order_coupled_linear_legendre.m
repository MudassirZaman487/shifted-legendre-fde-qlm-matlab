function [Xsol,Ysol] = solve_second_order_coupled_linear_legendre( ...
    X0,X1,Y0,Y1,F1,F2,a11,a12,a21,a22,sigma1,sigma2,eta,N,x)

B = shifted_legendre_fun_vector(N,eta,x);

Int1 = double(shifted_legendre_int(N,sigma1,eta));
Int2 = double(shifted_legendre_int(N,sigma2,eta));

h1 = X0 + X1*x;
h2 = Y0 + Y1*x;

H1 = double(shifted_legendre_coef(h1,N,eta,x));
H2 = double(shifted_legendre_coef(h2,N,eta,x));

F1 = double(shifted_legendre_coef(F1,N,eta,x));
F2 = double(shifted_legendre_coef(F2,N,eta,x));

V11 = double(shifted_legendre_var_matrix(a11,N,eta,x));
V12 = double(shifted_legendre_var_matrix(a12,N,eta,x));
V21 = double(shifted_legendre_var_matrix(a21,N,eta,x));
V22 = double(shifted_legendre_var_matrix(a22,N,eta,x));

B11 = eye(N) - Int1*V11;
B21 = -Int2*V12;

B12 = -Int1*V21;
B22 = eye(N) - Int2*V22;

SystemMat = [B11 B12; B21 B22];

RHS1 = H1*V11 + H2*V12 + F1;
RHS2 = H1*V21 + H2*V22 + F2;

RightSide = [RHS1 RHS2];

C = RightSide/SystemMat;

C1 = C(1:N);
C2 = C(N+1:2*N);

Xsol = vpa((C1*Int1 + H1)*B,8);
Ysol = vpa((C2*Int2 + H2)*B,8);

end
