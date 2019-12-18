function[Mi, y] = fc_ariane1(m)

global m_u
global k
global v
global delta_V

Mf = m_u + k(3)*m(3);%Mf,3
Mi = Mf + m(3);%Mi,3

y = v(3) * log(Mi/Mf);

Mf = Mi + k(2)*m(2);%Mf,2
Mi = Mf + m(2);%Mi,2

y = y + v(2) * log(Mi/Mf);

Mf = Mi + k(1)*m(1);%Mf,1
Mi = Mf + m(1);%Mi,1

y = y + v(1) * log(Mi/Mf);

y = y - delta_V;


end