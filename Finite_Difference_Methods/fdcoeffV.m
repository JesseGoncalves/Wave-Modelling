function c = fdcoeffV(k,xbar,x)
n = size(x,1);
A = ones(n,n);
xrow = (x(:)-xbar)';
for i = 2:n
    A(i,:) = (xrow.^(i-1))./factorial(i-1);
end
b = zeros(n,1);
b(k+1) = 1;
c = A\b;
c = c';
