function H_w_mat=MDOF_FRF_Point_Visc(M_mat,C_mat,K_mat,w)

H_w_mat=inv(-w^2*M_mat+1i*w*C_mat+K_mat);
