function [] = data()
%here we initialise some values for the whole problem

%radii, altitude
global R_t H_c R_c
R_t = 6378137;
H_c = 200e3;
R_c = R_t + H_c;

%aimed velocity
global mu V_c
mu = 3.986e14;
V_c = sqrt(mu/R_c);

%propelling velocity
global V_p
V_p = V_c;

%masses and properties of the propellants
global m_u m_e M_i k_ v_e alpha t_c
m_u = 1500;
m_e = [nan ; nan ; nan];
M_i = [nan ; nan ; nan ; m_u];
k_ = [0.1 ; 0.15 ; 0.2];
v_e = [2600 ; 3000 ; 4400];
alpha = [15 ; 10 ; 10];
t_c = [0 ; 0 ; 0];

%see the analytical part
global omega
for j = [1 2 3]
    omega(j) = k_(j)/(1 + k_(j));
end

end