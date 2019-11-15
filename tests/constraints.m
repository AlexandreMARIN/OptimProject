function y = constraints(x)

A = [1 0 ; 1 1 ; 0 1];

b = [0 ; 1 ; 0];

y = A*x+b;

end