function [H_cols,H_cols_SDOF_pages]=MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,w_col,m_vec,n_vec)

[P,Q]=size(EigVectors_Normalized);
if Q~=2*P,error('EigVectors_Normalized should be P x 2N matrix.'),end
N_w=size(w_col,1);
if ~isvector(m_vec),error('m_vec must be a vector!'),end
if (any(size(m_vec)~=size(n_vec)));error('Dimensions of m_row and n_row must be identical');end
N_FRF=length(m_vec);

H_cols_SDOF_pages=zeros(N_w,N_FRF,P);
for q=1:Q
    A_q_row=EigVectors_Normalized(m_vec,q).'.*EigVectors_Normalized(n_vec,q).';
    p=ceil(q/2);
    H_cols_SDOF_pages(:,:,p)=H_cols_SDOF_pages(:,:,p)+(1./(1i*w_col-EigValues_vec(q)))*A_q_row;
end
H_cols=sum(H_cols_SDOF_pages,3);