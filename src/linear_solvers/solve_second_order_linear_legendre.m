function approx_sol = solve_second_order_linear_legendre(y0,y1,G2,a0,b0,sigma,rho,eta,N,x)

B = shifted_legendre_fun_vector(N,eta,x);

IntMat = double(shifted_legendre_int(N,sigma,eta));
DiffMat = double(shifted_legendre_diff(N,rho,eta));

h = y0 + y1*x;

H = double(shifted_legendre_coef(h,N,eta,x));
F = double(shifted_legendre_coef(G2,N,eta,x));

VarA = double(shifted_legendre_var_matrix(a0,N,eta,x));
VarB = double(shifted_legendre_var_matrix(b0,N,eta,x));

OpMat = VarA + DiffMat*VarB;

SystemMat = eye(N) - IntMat*OpMat;
RightSide = H*OpMat + F;

C = RightSide/SystemMat;

approx_sol = vpa((C*IntMat + H)*B,8);

end
