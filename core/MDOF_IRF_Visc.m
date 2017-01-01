function h_cols= ...
MDOF_IRF_Visc(EigValues_mat, EigVectors_Normalized, t_column, n_row, m_row)

N=size(EigVectors_Normalized,1);
n=size(t_column,1);
n_col=size(n_row,2);

h_cols=zeros(n,n_col);
A_ind_row=sub2ind([N,N],n_row,m_row);

for r=1:2*N
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    A_r_temp_row=A_r(A_ind_row);
    h_cols=h_cols+exp(EigValues_mat(r,r)*t_column)*A_r_temp_row;
    
% %For Display only
% if imag(EigValues_mat(r,r))~=0 && mod(r,2)~=0
%     disp('IRF Response parameters');
%     r
%     2*abs(A_r)
% end
end

ind=find(abs(imag(h_cols))<=10*eps);
h_cols(ind)=real(h_cols(ind));