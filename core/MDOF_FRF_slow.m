function H_w_n_m_cols=MDOF_FRF_slow(MDOF_FRF_Point_func, w_column, N, n_row, m_row)

if any(size(n_row)~=size(m_row));error('Dimensions of n_row and m_row must be identical');end
i_col=size(n_row,2);

if (any((n_row>N)) || any((m_row>N))),error('Subscripts in n_row or m_row cannot exceed the system order N'),end

n_f=size(w_column,1);

H_ind_row=sub2ind([N,N],n_row,m_row);
H_w_n_m_cols=zeros(n_f,i_col);
for ii=1:n_f
    H_w_mat=MDOF_FRF_Point_func(w_column(ii));

    H_w_n_m_cols(ii,:)=H_w_mat(H_ind_row);
end
