function x_mul_m_vec= ...
SDOF_Harmonic_Response_Visc_mul_m(F0, w_0, w_n, zeta, t_vec, ...
								  ignoreTransient)    %Optional arguments

if nargin<6
    ignoreTransient=false;
end
r=w_0/w_n;
w_d=w_n*sqrt(1-zeta^2);

if (zeta ==0) && (r==1) && ~ignoreTransient
    x_mul_m_vec=-F0*w_n*t_vec.*cos(w_n*t_vec)/2/w_n^2;
else
    x_mul_m_vec=F0/w_n^2/((1-r^2)^2+(2*zeta*r)^2)*(-2*zeta*r*cos(w_0*t_vec)+(1-r^2)*sin(w_0*t_vec));
    if ~ignoreTransient
        x_mul_m_vec=x_mul_m_vec+F0/w_n^2/((1-r^2)^2+(2*zeta*r)^2)*(r*exp(-zeta*w_n*t_vec).*(2*zeta*cos(w_d*t_vec)+(2*zeta^2-1+r^2)*w_n/w_d*sin(w_d*t_vec)));
    end
end