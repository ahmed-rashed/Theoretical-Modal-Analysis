function x_rows= ...
MDOF_Free_Response_Visc(M,C, EigValues_mat, EigVectors_Normalized, x_0_col, x_dot_0_col, t_row)

N=size(EigVectors_Normalized,1);
n_col=size(t_row,2);
if size(x_0_col,1)~=N,error('initial conditions must have the same length as the order of the square matrix M '),end

x_rows=zeros(N,n_col);
col1=M*x_0_col;
col2=C*x_0_col+M*x_dot_0_col;
for r=1:2*N
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    x_rows=x_rows+EigValues_mat(r,r)*A_r*col1*exp(EigValues_mat(r,r)*t_row)+A_r*col2*exp(EigValues_mat(r,r)*t_row);
    
%For Display only
if imag(EigValues_mat(r,r))~=0 && mod(r,2)~=0
    disp('Free Response parameters');
    r
    2*(real(EigValues_mat(r,r))*real(A_r)-abs(imag(EigValues_mat(r,r)))*imag(A_r))*col1+2*real(A_r)*col2
    -2*(abs(imag(EigValues_mat(r,r)))*real(A_r)+real(EigValues_mat(r,r))*imag(A_r))*col1-2*imag(A_r)*col2
end

end

x_rows(imag(x_rows)<100*eps)=real(x_rows(imag(x_rows)<100*eps));