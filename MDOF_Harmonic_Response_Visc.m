function x_rows=MDOF_Harmonic_Response_Visc(EigValues_mat, EigVectors_Normalized, F0_col, w_0_col, t_row,ignoreTransient)

N=size(EigVectors_Normalized,1);
n_col=size(t_row,2);
n_row=size(F0_col,1);

[w_r_vec, zeta_r_vec]=MDOF_Modal_Param_Visc(EigValues_mat);
if any(isnan(w_r_vec))
    warning('Real pole encountered. This function assumes all poles are complex conjugates.');
end

x_rows=zeros(n_row,n_col);
semiCol1=zeros(N,n_col);
semiCol2=zeros(N,n_col);
for r=1:2:2*N-1
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    for jj=1:N
        semiCol1(jj,:)=SDOF_Harmonic_Response_dot_Visc_mul_m(F0_col(jj), w_0_col(jj), w_r_vec((r+1)/2), zeta_r_vec((r+1)/2), t_row,ignoreTransient);
        semiCol2(jj,:)=SDOF_Harmonic_Response_Visc_mul_m(F0_col(jj), w_0_col(jj), w_r_vec((r+1)/2), zeta_r_vec((r+1)/2), t_row,ignoreTransient);
    end
    x_rows=x_rows+real(A_r)*semiCol1-(real(EigValues_mat(r,r))*real(A_r)+imag(EigValues_mat(r,r))*imag(A_r))*semiCol2;
end

x_rows=2*x_rows;

x_rows(imag(x_rows)<100*eps)=real(x_rows(imag(x_rows)<100*eps));