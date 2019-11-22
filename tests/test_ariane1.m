addpath("../")
clearvars

%parameters of the problem
global v
global k
global m_u
global delta_V

v = [2647.2 ; 2922.4 ; 4344.3];
k = [0.1101 ; 0.1532 ; 0.2154];
m_u = 1700;
delta_V = 11527;

m_star = [145349 ; 31215 ; 7933];
f_star = 208691;

fprintf("f(m_star) = %f\n\nc(m_star) = %f\n\n", f_ariane1(m_star), c_ariane1(m_star))

m_inf = [50000 ; 10000 ; 1000];
m_sup = [500000 ; 50000 ; 10000];
m0 = (m_inf + m_sup)/2;
mindm = 1;
mindf = 1;
maxnbcall = 90;
maxnbiter = 150;
tol = 0.01;
[m, lambda, iter, normGradLag] = SQP(m0, @f_ariane1, @c_ariane1, tol, m_inf, m_sup, maxnbiter, maxnbcall, mindm, mindf);