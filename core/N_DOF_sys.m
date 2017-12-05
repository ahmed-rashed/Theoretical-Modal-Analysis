function [M_mat,C_mat,K_mat]=N_DOF_sys(m_vec,c_vec,k_vec)
N=length(m_vec);

M_mat=diag(m_vec);
C_mat=diag(c_vec(1:N)+c_vec(2:N+1))+diag(-c_vec(2:N),1)+diag(-c_vec(2:N),-1);
K_mat=diag(k_vec(1:N)+k_vec(2:N+1))+diag(-k_vec(2:N),1)+diag(-k_vec(2:N),-1);