function H_cols=MDOF_FRF_slow(MDOF_FRF_Point_func,w_column,N,n_row,m_row)

if any(size(n_row)~=size(m_row));error('Dimensions of n_row and m_row must be identical');end
N_cols=size(n_row,2);

if (any((n_row>N)) || any((m_row>N))),error('Subscripts in n_row or m_row cannot exceed the system order N'),end

N_w=size(w_column,1);

H_ind_row=sub2ind([N,N],n_row,m_row);
H_cols=nan(N_w,N_cols);
for n_w=1:N_w
    H_mat=MDOF_FRF_Point_func(w_column(n_w));

    H_cols(n_w,:)=H_mat(H_ind_row);
end
