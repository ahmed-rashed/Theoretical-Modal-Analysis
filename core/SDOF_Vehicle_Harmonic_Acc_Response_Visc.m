function y_Acc_vec=SDOF_Vehicle_Harmonic_Acc_Response_Visc ...
						 (Y_0,w_0,w_n,zeta,t_vec, ...
                                              ignoreTransient)    %Optional arguments

if nargin<6
    ignoreTransient=false;
end

r_0=w_0/w_n;
if (zeta ==0) && (r_0==1) && ~ignoreTransient
    y_Acc_vec=Y_0*w_n^3/2*t_vec.*cos(w_n*t_vec);
else
    y_Acc_vec=2*zeta*r_0^5*cos(w_0*t_vec)-r_0^2*((2*zeta*r_0)^2+1-r_0^2)*sin(w_0*t_vec);
    if ~ignoreTransient
        w_d=w_n*sqrt(1-zeta^2);
        y_Acc_vec=y_Acc_vec+exp(-zeta*w_n*t_vec).*(-2*zeta*r_0^3*(1+zeta/w_n-zeta^2)*cos(w_d*t_vec)+r_0*(1-r_0^2*(1+2*zeta^2))*(w_n*(1-zeta^2)+zeta)*sin(w_d*t_vec)/w_d);
    end

    y_Acc_vec=Y_0*w_n^2/((1-r_0^2)^2+(2*zeta*r_0)^2)*y_Acc_vec;
end