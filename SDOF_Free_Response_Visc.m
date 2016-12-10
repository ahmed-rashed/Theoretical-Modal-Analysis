function x_vec=SDOF_Free_Response_Visc(w_n, zeta, x0, x_dot_0, t_vec)

if zeta~=1
    w_d=w_n*sqrt(1-zeta^2);
    x_vec=exp(-zeta*w_n*t_vec).*(x0*cos(w_d*t_vec)+(zeta*w_n*x0+x_dot_0)*sin(w_d*t_vec)/w_d);
else
    x_vec=exp(-w_n*t_vec).*(x0+(w_n*x0+x_dot_0)*t_vec);
end