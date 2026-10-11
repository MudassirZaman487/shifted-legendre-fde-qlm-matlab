function VarMat = shifted_legendre_var_matrix(g,N,eta,x)

C = shifted_legendre_coef(g,N,eta,x);

for l = 0:N-1
    for n = 0:N-1
        for k = 0:N-1
            P(k+1) = ((2*n+1)/eta)*C(k+1)* ...
                     three_shifted_legendre_product(k,l,n,eta);
        end
        D(l+1,n+1) = sum(P);
    end
end

VarMat = double(D);

end
