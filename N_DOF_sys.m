function [M,C,K]=N_DOF_sys(m_vec,c_vec,k_vec)
N=length(m_vec);

M=diag(m_vec);
C=diag(c_vec(1:N)+c_vec(2:N+1))+diag(-c_vec(2:N),1)+diag(-c_vec(2:N),-1);
K=diag(k_vec(1:N)+k_vec(2:N+1))+diag(-k_vec(2:N),1)+diag(-k_vec(2:N),-1);