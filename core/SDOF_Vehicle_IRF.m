function h_vec=SDOF_Vehicle_IRF(w_n,zeta,t_vec)

if zeta==1
    h_vec=w_n*exp(-w_n*t_vec).*(2-w_n*t_vec);
else
    w_d=w_n*sqrt(1-zeta^2);
    h_vec=w_n*exp(-zeta*w_n*t_vec).*(2*zeta*cos(w_d*t_vec)+(1-2*zeta^2)/sqrt(1-zeta^2)*sin(w_d*t_vec));
end

h_vec(t_vec<0)=0;