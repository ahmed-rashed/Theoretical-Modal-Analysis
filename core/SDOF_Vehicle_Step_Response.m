function h_vec=SDOF_Vehicle_Step_Response(Y_0,w_n, zeta, t_vec)

if zeta==1
    h_vec=1+exp(-w_n*t_vec).*(w_n*t_vec-1);
else
    w_d=w_n*sqrt(1-zeta^2);
    h_vec=1+exp(-zeta*w_n*t_vec).*(zeta*w_n*sin(w_d*t_vec)/w_d-cos(w_d*t_vec));
end
h_vec(t_vec<0)=0;
h_vec=Y_0*h_vec;