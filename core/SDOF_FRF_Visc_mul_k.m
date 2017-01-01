function H_k_vec= SDOF_FRF_Visc_mul_k(r_vec,zeta)

H_k_vec=1./(1-r_vec.^2+2*1i*zeta*r_vec);