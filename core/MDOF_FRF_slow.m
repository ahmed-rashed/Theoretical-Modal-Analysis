function H_cols=MDOF_FRF_slow(MDOF_FRF_Point_func,w_column,N,m_row,n_row)

if any(size(m_row)~=size(n_row));error('Dimensions of m_row and n_row must be identical');end
N_cols=size(m_row,2);

if (any((m_row>N)) || any((n_row>N))),error('Subscripts in m_row or n_row cannot exceed the system order N'),end

N_w=size(w_column,1);

H_ind_row=sub2ind([N,N],m_row,n_row);
H_cols=nan(N_w,N_cols);
for n_w=1:N_w
    H_mat=MDOF_FRF_Point_func(w_column(n_w));

    H_cols(n_w,:)=H_mat(H_ind_row);
end
