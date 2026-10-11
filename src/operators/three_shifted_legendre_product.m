function T = three_shifted_legendre_product(a,b,c,eta)
syms l m n
A = @(i,k) ((-1)^(i+k)*gamma(i+k+1))/( ...
           gamma(i-k+1)*(eta^k)*(gamma(k+1)^2));
QQ = 0;

for l = 0:a
    for m = 0:b
        for n = 0:c
            aa = A(a,l);
            bb = A(b,m);
            cc = A(c,n);
            integral_part = eta^(l+m+n+1)/(l+m+n+1);
            QQ = QQ + aa*bb*cc*integral_part;
        end
    end
end

T = double(QQ);
