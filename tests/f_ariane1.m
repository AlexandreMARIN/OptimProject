function[M_i_1] = f_ariane1(m)

global m_u
global k

M_i_1 = m(1)*(1+k(1)) + m(2)*(1+k(2)) + m(3)*(1+k(3)) + m_u;
end