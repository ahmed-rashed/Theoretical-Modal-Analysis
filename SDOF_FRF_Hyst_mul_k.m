function H_k_vec= SDOF_FRF_Hyst_mul_k(r_vec,eta)

H_k_vec=1./(1-r_vec.^2+i*eta);