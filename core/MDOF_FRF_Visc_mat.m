function H_cols_pages=MDOF_FRF_Visc_mat(s_q_vec,A_q_mat_pages,w_col)

[M,N,Q]=size(A_q_mat_pages);
if length(s_q_vec)~=Q,error('3rd dimension of A_q_mat_pages must be the same as the length of s_q_vec!'),end
N_w=size(w_col,1);

H_cols_pages=zeros(N_w,M,N);
for q=1:Q
    H_cols_pages=H_cols_pages+permute(A_q_mat_pages(:,:,q),[3,1,2])./(1i*w_col-s_q_vec(q));
end