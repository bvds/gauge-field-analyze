% EXAMPLE 3: A is diagonal, singular and indefinite.
h = 1;  a = -10; b = -a; n = 2*b/h + 1;
A = spdiags((a:h:b)', 0, n, n);
b = ones(n,1);  rtol   = 1e-6;    shift = 0;    maxxnorm = 1e2;
M = [];         Acondlim = [];    show  = true;
x = minresqlp(A,b,rtol,n,M,shift,maxxnorm,Acondlim,show);
