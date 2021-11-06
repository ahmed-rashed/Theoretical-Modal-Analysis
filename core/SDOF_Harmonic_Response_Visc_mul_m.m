function x_mul_m_vec= ...
SDOF_Harmonic_Response_Visc_mul_m(F0,w_0,w_n,zeta,t_vec, ...
								  ignoreTransient)    %Optional arguments

if nargin<6
    ignoreTransient=false;
end
r_0=w_0/w_n;

if (zeta ==0) && (r_0==1) && ~ignoreTransient
    x_mul_m_vec=-F0/2/w_n*t_vec.*cos(w_n*t_vec);
else
    x_mul_m_vec=-2*zeta*r_0*cos(w_0*t_vec)+(1-r_0^2)*sin(w_0*t_vec);
    if ~ignoreTransient
        w_d=w_n*sqrt(1-zeta^2);
        x_mul_m_vec=x_mul_m_vec+r_0*exp(-zeta*w_n*t_vec).*(2*zeta*cos(w_d*t_vec)+(2*zeta^2-1+r_0^2)*w_n/w_d*sin(w_d*t_vec));
    end

    x_mul_m_vec=F0/w_n^2/((1-r_0^2)^2+(2*zeta*r_0)^2)*x_mul_m_vec;
end