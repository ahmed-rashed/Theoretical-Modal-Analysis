function [H_cols,H_cols_SDOF_pages]=MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,w_col,m_vec,n_vec)

[N,Q]=size(EigVectors_Normalized);
if Q~=2*N,error('EigVectors_Normalized should be N x 2N matrix.'),end
N_w=size(w_col,1);
if ~isvector(m_vec),error('m_vec must be a vector!'),end
if (any(size(m_vec)~=size(n_vec)));error('Dimensions of m_row and n_row must be identical');end
N_FRF=length(m_vec);

H_cols_SDOF_pages=zeros(N_w,N_FRF,N);
for q=1:2*N
    A_r_row=EigVectors_Normalized(m_vec,q).'.*EigVectors_Normalized(n_vec,q).';
    r=ceil(q/2);
    H_cols_SDOF_pages(:,:,r)=H_cols_SDOF_pages(:,:,r)+(1./(1i*w_col-EigValues_vec(q)))*A_r_row;
end
H_cols=sum(H_cols_SDOF_pages,3);