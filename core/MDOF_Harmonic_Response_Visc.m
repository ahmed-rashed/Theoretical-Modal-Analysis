function x_rows=MDOF_Harmonic_Response_Visc(EigValues_vec,EigVectors_Normalized,F0_col,w_0_col,t_row,ignoreTransient)

P=size(EigVectors_Normalized,1);
n_col=size(t_row,2);
n_row=size(F0_col,1);

[w_p_vec,zeta_p_vec]=pole2modal_visc(EigValues_vec);
if any(isnan(w_p_vec))
    warning('Real pole encountered. This function assumes all poles are complex conjugates.');
end

x_rows=zeros(n_row,n_col);
semiCol1=zeros(P,n_col);
semiCol2=zeros(P,n_col);
for p=1:2:2*P-1
    A_p=EigVectors_Normalized(:,p)*EigVectors_Normalized(:,p).';
    for jj=1:P
        semiCol1(jj,:)=SDOF_Harmonic_Response_dot_Visc_mul_m(F0_col(jj),w_0_col(jj),w_p_vec((p+1)/2),zeta_p_vec((p+1)/2),t_row,ignoreTransient);
        semiCol2(jj,:)=SDOF_Harmonic_Response_Visc_mul_m(F0_col(jj),w_0_col(jj),w_p_vec((p+1)/2),zeta_p_vec((p+1)/2),t_row,ignoreTransient);
    end
    x_rows=x_rows+real(A_p)*semiCol1-(real(EigValues_vec(p))*real(A_p)+imag(EigValues_vec(p))*imag(A_p))*semiCol2;
end

x_rows=2*x_rows;

x_rows(imag(x_rows)<100*eps)=real(x_rows(imag(x_rows)<100*eps));