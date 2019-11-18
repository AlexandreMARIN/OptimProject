function [x, lambda] = SQP(x0, f, c, epsilon, nbiter)
%the SQP algorithm
%Input:
%x0 : first estimate of a minimiser
%f, c : objective, constraints
%epsilon : tolerance
%nbiter : maximal number of iterations
%Output:
%x : an estimate of a minimiser
%lambda : multipliers of Lagrange

[n, ~] = size(x0);
h = 1e-3;%for computing gradients
c1 = 0.1;%see Armijo
nbit_Armijo = 100;%number of iterations for Armijo

k = 0;%rank of the current iteration
x = x0;



%compute A, Q, g, b
H = eye(n);
fx = f(x);
cx = c(x);
[g, A] = gradient(x, fx, cx, h, f, c);
b = -cx;

%we solve the quadratic problem
[d_QP, lambda] = KKT_quad_lin(H, g, A, b);
new_x = x + d_QP;

%globalisation
Fx = merit(fx, cx, rho);
prev_x = x;%previous value of x
rho = max(lambda)+1;
[x, fx, cx, Fx] = Armijo(new_x, x, Fx, g, cx, c1, f, c, rho, nbit_Armijo);


while k < nbiter && Fx > epsilon
    
    k = k + 1;
    
    %compute A, Q, g, b
    [g, A] = gradient(x, fx, cx, h, f, c);
    b = -cx;
    H = compute_H(H, prev_x, x, lambda, f, c);

    %we solve the quadratic problem
    [d_QP, lambda] = KKT_quad_lin(H, g, A, b);
    new_x = x + d_QP;
    
    %globalisation
    prev_x = x;
    rho = max(lambda)+1;
    [x, fx, cx, Fx] = Armijo(new_x, x, Fx, g, cx, c1, f, c, rho, nbit_Armijo);

end

end