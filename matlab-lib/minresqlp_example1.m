%clear functions
% minresqlp

n = 10;                                  e = ones(n,1);
A = spdiags([-2*e 4*e -2*e],-1:1,n,n);    M = spdiags(4*e,0,n,n);
b = sum(A,2);             rtol = 1e-10;   maxit = 50;

% Alternatively, use a function for A.
%A = @(x)Afun(x,n);

x = minresqlp(A,b,rtol,maxit,M);
x

function y = Afun(x,n)
    y = 4*x;
    y(2:n)   = y(2:n)   - 2*x(1:n-1);
    y(1:n-1) = y(1:n-1) - 2*x(2:n);
end
