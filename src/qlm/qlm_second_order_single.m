function U = qlm_second_order_single(N,sigma,rho,eta,y0,y1,f,a_x,b_x,ITT,x)

U = cell(ITT,1);

U{1} = solve_second_order_linear_legendre(y0,y1,f,a_x,b_x,sigma,rho,eta,N,x);
U{1} = expand(vpa(U{1},8));

for i = 1:ITT-1

    Ui = expand(vpa(U{i},8));

    an = expand(a_x + 2*Ui);
    bn = b_x;
    fn = expand(f - Ui^2);

    U{i+1} = solve_second_order_linear_legendre(y0,y1,fn,an,bn,sigma,rho,eta,N,x);
    U{i+1} = expand(vpa(U{i+1},8));

    disp(['QLM iteration = ',num2str(i+1),', N = ',num2str(N), ...
          ', sigma = ',num2str(sigma),', rho = ',num2str(rho)])

end

end
