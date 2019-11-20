function [x, lambda,k] = SQP(x0, f, c, epsilon, nbiter, x_inf, x_sup)
%the SQP algorithm
%Input:
%x0 : first estimate of a minimiser
%f, c : objective, constraints
%epsilon : tolerance
%nbiter : maximal number of iterations
%Output:
%x : an estimate of a minimiser
%lambda : multipliers of Lagrange
x = x0;

[n, ~] = size(x);
h = 1e-3;%for computing gradients
c1 = 0.1;%see Armijo
nbit_Armijo = 1000;%number of iterations for Armijo



k = 0;%rank of the current iteration




%compute A, Q, g, b
H = eye(n);
fx = f(x);
cx = c(x);
[g, A] = gradient(x, fx, cx, h, f, c);
b = -cx;

%we solve the quadratic problem
[d_QP, lambda] = KKT_quad_lin(H, g, A, b);
fprintf("LAMBDA : \n");
disp(lambda);
fprintf("END LAMBDA\n")
new_x = x + d_QP;
ind_inf = new_x < x_inf;
new_x(ind_inf) = x_inf(ind_inf);
ind_sup = new_x > x_sup;
new_x(ind_sup) = x_sup(ind_sup);


%globalisation

rho = max(lambda)*100;
fprintf("RHO : \n");
disp(rho);
fprintf("END RHO ");


Fx = merit(fx, cx, rho);
prev_x = x;%previous value of x
fprev = fx;
cprev = cx;
[x, fx, cx, Fx] = Armijo(new_x, x, Fx, g, cx, c1, f, c, rho, nbit_Armijo);


while (k < nbiter) && (Fx > epsilon)
    
    k = k + 1;
    
    %compute A, Q, g, b
    [g, A] = gradient(x, fx, cx, h, f, c);
    b = -cx;
    H = compute_H(H, prev_x, fprev, cprev, x, fx, cx, lambda, f, c);

    %we solve the quadratic problem
    [d_QP, lambda] = KKT_quad_lin(H, g, A, b);
    new_x = x + d_QP;
    ind_inf = new_x < x_inf;
    new_x(ind_inf) = x_inf(ind_inf);
    ind_sup = new_x > x_sup;
    new_x(ind_sup) = x_sup(ind_sup);
    
    %globalisation
    prev_x = x;
    fprev = fx;
    cprev = cx;
    rho = max(lambda)*100;
    [x, fx, cx, Fx] = Armijo(new_x, x, Fx, g, cx, c1, f, c, rho, nbit_Armijo);

end

end