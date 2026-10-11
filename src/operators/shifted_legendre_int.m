function IntMat = shifted_legendre_int(N,zeta,eta)

syms k l

for i = 0:N-1
    for j = 0:N-1
        aa = ((-1)^(i+k)*gamma(i+k+1))/( ...
             gamma(i-k+1)*(eta^k)*(gamma(k+1)^2));
        bb = gamma(k+1)/gamma(k+zeta+1);
        cc = ((2*j+1)/eta)* ...
             (((-1)^(j+l)*gamma(j+l+1))/( ...
             gamma(j-l+1)*(eta^l)*(gamma(l+1)^2))) * ...
             (eta^(k+zeta+l+1)/(k+zeta+l+1));
        P(i+1,j+1) = symsum(aa*bb*symsum(cc,l,0,j),k,0,i);
    end
end
IntMat = double(P);
end
