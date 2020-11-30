function [H_cols,H_cols_SDOF_layers]=MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,w_col,n_row,m_row)

[N,Q]=size(EigVectors_Normalized);
if Q~=2*N,error('EigVectors_Normalized should be N x 2N matrix.'),end
N_w=size(w_col,1);
if (any(size(n_row)~=size(m_row)));error('Dimensions of n_row and m_row must be identical');end
N_FRF=size(n_row,2);

H_cols_SDOF_layers=zeros(N_w,N_FRF,N);
for q=1:2*N
    A_r_row=EigVectors_Normalized(n_row,q).'.*EigVectors_Normalized(m_row,q).';
    r=ceil(q/2);
    H_cols_SDOF_layers(:,:,r)=H_cols_SDOF_layers(:,:,r)+(1./(1i*w_col-EigValues_vec(q)))*A_r_row;
end
H_cols=sum(H_cols_SDOF_layers,3);