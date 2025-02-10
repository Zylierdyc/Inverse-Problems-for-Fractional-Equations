function g = coe_a(tau,alpha,sigma,N)

% 求积系数
a(1) = sigma^(1-alpha);
for i = 1:N-1
    a(i+1) = (i+sigma)^(1-alpha)-(i-1+sigma)^(1-alpha);
    b(i+1) = 1/(2-alpha)*((i+sigma)^(2-alpha)-(i-1+sigma)^(2-alpha)) ...
        -1/2*((i+sigma)^(1-alpha)+(i-1+sigma)^(1-alpha));
end
for  i = 2:N-1
    g(i) = tau^(-alpha)/gamma(2-alpha)*(a(N-i+1)+b(N-i+2)-b(N-i+1));
end
if N == 1
    g(N) = tau^(-alpha)/gamma(2-alpha)*a(1);
else
    g(N) = tau^(-alpha)/gamma(2-alpha)*(a(1)+b(2));
    g(1) = tau^(-alpha)/gamma(2-alpha)*(a(N)-b(N)); 
end
end