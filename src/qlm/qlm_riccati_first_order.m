%% ==========================================================
%  QLM solver for Riccati equation
%
%  D^gamma u = -u^2 + 1
%  u(0) = u0
%
%  Linearization:
%  -u_{n+1}^2 is replaced by
%  -2*u_n*u_{n+1} + u_n^2
%
%  Therefore:
%  D^gamma u_{n+1} = a_n(x)*u_{n+1} + r_n(x)
%
%  where
%  a_n(x) = -2*u_n(x)
%  r_n(x) = u_n^2(x) + 1
% ==========================================================

function U = qlm_riccati_first_order(N,gamma,eta,u0,ITT,x)

U = cell(ITT,1);

a_start = sym(0);
r_start = sym(1);

U{1} = solve_first_order_linear_legendre(u0,r_start,a_start,gamma,eta,N,x);
U{1} = expand(vpa(U{1},8));

for i = 1:ITT-1

    Ui = expand(vpa(U{i},8));

    a_i = expand(-2*Ui);
    r_i = expand(Ui^2 + 1);

    U{i+1} = solve_first_order_linear_legendre(u0,r_i,a_i,gamma,eta,N,x);
    U{i+1} = expand(vpa(U{i+1},8));

    disp(['QLM iteration = ',num2str(i+1),', N = ',num2str(N), ...
          ', gamma = ',num2str(gamma)])

end

end
