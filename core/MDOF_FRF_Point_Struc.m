function H_w_mat=MDOF_FRF_Point_Struc(M_mat, D_mat, K_mat, w)

H_w_mat=inv(K_mat+1i*D_mat-w^2*M_mat);
