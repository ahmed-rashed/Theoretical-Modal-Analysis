function x_mul_m_dot_vec=SDOF_Harmonic_Response_Visc_mul_m_dot ...
						 (F0, w_0, w_n, zeta, t_vec, ...
                                              ignoreTransient)    %Optional arguments

if nargin<6
    ignoreTransient=false;
end

r=w_0/w_n;
w_d=w_n*sqrt(1-zeta^2);
if (zeta ==0) && (r==1) && ~ignoreTransient
    x_mul_m_dot_vec=F0*w_0*t_vec/(w_0+w_n).*sin((w_0+w_n)*t_vec/2);
else
    x_mul_m_dot_vec=F0*r/w_n*((1-r^2)*cos(w_0*t_vec)+2*zeta*r*sin(w_0*t_vec))/((1-r^2)^2+(2*zeta*r)^2);
    if ~ignoreTransient
        x_mul_m_dot_vec=x_mul_m_dot_vec-F0*r/w_n*exp(-zeta*w_n*t_vec).*((1-r^2)*cos(w_d*t_vec)+zeta*w_n*(r^2+1)*sin(w_d*t_vec)/w_d)/((1-r^2)^2+(2*zeta*r)^2);
    end
end