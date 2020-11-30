function x_mul_m_dot_vec=SDOF_Harmonic_Response_dot_Visc_mul_m ...
						 (F0, w_0, w_n, zeta, t_vec, ...
                                              ignoreTransient)    %Optional arguments

if nargin<6
    ignoreTransient=false;
end

r_0=w_0/w_n;
w_d=w_n*sqrt(1-zeta^2);
if (zeta ==0) && (r_0==1) && ~ignoreTransient
    x_mul_m_dot_vec=F0*w_0*t_vec/(w_0+w_n).*sin((w_0+w_n)*t_vec/2);
else
    x_mul_m_dot_vec=F0*r_0/w_n*((1-r_0^2)*cos(w_0*t_vec)+2*zeta*r_0*sin(w_0*t_vec))/((1-r_0^2)^2+(2*zeta*r_0)^2);
    if ~ignoreTransient
        x_mul_m_dot_vec=x_mul_m_dot_vec-F0*r_0/w_n*exp(-zeta*w_n*t_vec).*((1-r_0^2)*cos(w_d*t_vec)+zeta*w_n*(r_0^2+1)*sin(w_d*t_vec)/w_d)/((1-r_0^2)^2+(2*zeta*r_0)^2);
    end
end