function[Mi, y] = stageproblem(m)
%it gives the objective and the constaints for the stage problem

global m_u k_ v_e V_p

Mf = m_u + k_(3)*m(3);%Mf,3
Mi = Mf + m(3);%Mi,3

y = v_e(3) * log(Mi/Mf);

Mf = Mi + k_(2)*m(2);%Mf,2
Mi = Mf + m(2);%Mi,2

y = y + v_e(2) * log(Mi/Mf);

Mf = Mi + k_(1)*m(1);%Mf,1
Mi = Mf + m(1);%Mi,1

y = y + v_e(1) * log(Mi/Mf);

y = y - V_p;


end