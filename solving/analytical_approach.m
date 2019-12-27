%we solve analytically the stage problem, see the questions 1.4.1-1.4.4
%we use a Newton algorithm, the function definitions at the end of that file
%match the functions we use theoretically
clearvars

global m_u M_i m_e v_e k_ omega V_p

%we load some variables
data()

%remember that in data(), there is this assignment: V_p = V_c
V_p = V_p*1.2;

x30 = 3;%good initial value, see the plotted graph
tol = 1e-8;%tolerance
maxnbiter = 10;%maximum number of iterations

%Newton algorithm
x3 = x30;
for iter = 1:maxnbiter
    cx = c(x3);
    if abs(cx) < tol
        break
    end
    x3 = x3 - cx/cprime(x3);
end

%we plot a graph for you to understand the choice of the initial value
val = linspace(2.6, 3.5, 100);
cval = zeros(100,1);
for i = 1:100
    cval(i) = c(val(i));
end
plot(val, cval, '-r', x3, cx, '+b')
legend(["constraint", "estimate"])
title("check the estimate of $x_{3}$", 'Interpreter', 'latex')

%we construct x and the masses of our problem
x = [0 0 x3];
for j = [2 1]
    x(j) = (1 - (v_e(3)/v_e(j))*(1 - omega(3)*x(3)))/omega(j);
end
for j = [3 2 1]
    M_i(j) = M_i(j+1)/((1+k_(j))/x(j) - k_(j));
    m_e(j) = (M_i(j) - M_i(j+1))/(1+k_(j));
end

%this is the asked multiplier of Lagrange
lambda = -m_u/(M_i(1)*v_e(3)*(1 - omega(3)*x(3)));

%here we check that the formula given by the KKT theorem
%the other two equations are equivalent
KKT = ((1+k_(1))/x(1) - k_(1))*((1+k_(2))/x(2) - k_(2))*((1+k_(3))/(x(3)^2))...
    + lambda*v_e(3)/x(3);

fprintf("equation given by the KKT theorem : %f\n", KKT)

%function definitions
function [cx] = c(x)
global v_e V_p omega
cx = v_e(3)*log(x) + v_e(2)*log((1 - (v_e(3)/v_e(2))*(1 - omega(3)*x))/omega(2))...
   + v_e(1)*log((1 - (v_e(3)/v_e(1))*(1 - omega(3)*x))/omega(1)) - V_p;
end

function [cpx] = cprime(x)
global v_e omega
cpx = v_e(3)*(1/x + omega(3)/(1 - (v_e(3)/v_e(2))*(1 - omega(3)*x))...
    + omega(3)/(1 - (v_e(3)/v_e(1))*(1 - omega(3)*x)));
end