% EXAMPLE 4: Use this matrix-vector product function below
% Call MINRESQLP using the anonymous function handle:

b = rand(400,1);
x = minresqlp(@(x,N)Afun(x,20), b);

% NOTE: The first x is vector output. The 'x''s in the first input are
% abstraction for the input of fun. Note how specific values for additional
% input parameters of Afun (such as N) is passed into the function.
%
% Alternatively, we can define the function handle
%       g = @(x,N)Afun(x,20);
% and then call MINRESQLP using
%       x = minresqlp(g, b);
% or
%       x = minresqlp('@(x,N)Afun(x,20)', b);

function y = Afun(x, N)

n = length(x);

e = ones(N,1);
T = spdiags([e e e], -1:1, N, N);

M = n/N;
y = zeros(n,1);

for i=1:M
  i1 = (i-2)*N + 1;
  i2 = (i-1)*N + 1;
  i3 = i    *N + 1;
  i4 = i3 + N - 1;
  if (i1 >=1 && i3 <= n)
    y(i2:(i3-1)) = T * (x(i1:(i2-1)) + x(i2:(i3-1)) + x(i3:i4));
  elseif (i1 <= 0 && i3 <= n)
    y(i2:(i3-1)) = T * (               x(i2:(i3-1)) + x(i3:i4));
  elseif (i1 >=1 && i3 > n)
    y(i2:(i3-1)) = T * (x(i1:(i2-1)) + x(i2:(i3-1))            );
  end
end
end

