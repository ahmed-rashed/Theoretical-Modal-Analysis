function h_m_vec=SDOF_IRF_Visc_mul_m(w_n, zeta, t_vec)

if zeta~=1
    w_d=w_n*sqrt(1-zeta^2);
    h_m_vec=exp(-zeta*w_n*t_vec).*sin(w_d*t_vec)/w_d;
else
    h_m_vec=exp(-w_n*t_vec).*t_vec;
end