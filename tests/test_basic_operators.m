% BASIC OPERATOR SMOKE TEST (not yet executed in MATLAB).
% Run setup_paths from the repository root before using this script.
syms x
N = 4; eta = 1;
B = shifted_legendre_fun_vector(N, eta, x);
C = shifted_legendre_coef(x, N, eta, x);
I = shifted_legendre_int(N, 1, eta);
D = shifted_legendre_diff(N, 1, eta);
V = shifted_legendre_var_matrix(x, N, eta, x);
xcheck = 0.3;
intResult = double(subs(C*I*B, x, xcheck));
diffResult = double(subs(C*D*B, x, xcheck));
C1 = shifted_legendre_coef(sym(1), N, eta, x);
varResult = double(subs(C1*V*B, x, xcheck));
assert(abs(intResult - xcheck^2/2) < 1e-7, 'Integral test failed');
assert(abs(diffResult - 1) < 1e-7, 'Derivative test failed');
assert(abs(varResult - xcheck) < 1e-7, 'Product matrix test failed');
disp('Basic operator tests passed.');
