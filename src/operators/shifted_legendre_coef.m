function C = shifted_legendre_coef(f,N,eta,x)

A = shifted_legendre_fun_vector(N,eta,x);

F = int(f*A,x,0,eta);

for i = 0:N-1
    B(i+1) = double((2*i+1)/eta)*F(i+1);
end

C = B;

end
