% EXAMPLE 2: A is Laplacian on a 50 by 50 grid, singular and indefinite.
n = 50;  N = n^2;  e = ones(n,1);
B = spdiags([e e e], -1:1, n, n);
A = sparse([],[],[],N,N,(3*n-2)^2);
for i=1:n
    A((i-1)*n+1:i*n,(i-1)*n+1:i*n) = B;
    if i*n+1 < n*n,   A(i*n+1:(i+1)*n,(i-1)*n+1:i*n)     = B; end
    if (i-2)*n+1 > 0, A((i-2)*n+1:(i-1)*n,(i-1)*n+1:i*n) = B; end
end
b = sum(A,2);   rtol   = 1e-5;    shift = 0;    maxxnorm = 1e2;
M = [];         Acondlim = [];    show  = true;
x = minresqlp(A,b,rtol,N,M,shift,maxxnorm,Acondlim,show);
