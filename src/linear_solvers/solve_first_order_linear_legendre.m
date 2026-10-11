function approx_sol = solve_first_order_linear_legendre(y0,G2,a0,gamma,eta,N,x)

B = shifted_legendre_fun_vector(N,eta,x);

IntMat = double(shifted_legendre_int(N,gamma,eta));

H = double(shifted_legendre_coef(y0,N,eta,x));
F = double(shifted_legendre_coef(G2,N,eta,x));

VarMat = double(shifted_legendre_var_matrix(a0,N,eta,x));

SystemMat = eye(N) - IntMat*VarMat;
RightSide = H*VarMat + F;

C = RightSide/SystemMat;

approx_sol = vpa((C*IntMat + H)*B,8);

N

end
