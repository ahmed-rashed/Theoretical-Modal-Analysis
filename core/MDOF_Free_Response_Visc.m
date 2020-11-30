function x_rows= ...
MDOF_Free_Response_Visc(M,C,EigValues_vec,EigVectors_Normalized,x_0_col,x_dot_0_col,t_row)

N=size(EigVectors_Normalized,1);
n_col=size(t_row,2);
if size(x_0_col,1)~=N,error('initial conditions must have the same length as the order of the square matrix M '),end

x_rows=zeros(N,n_col);
v_1=M*x_0_col;
v_3=C*x_0_col+M*x_dot_0_col;
for r=1:2*N
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    x_rows=x_rows+EigValues_vec(r)*A_r*v_1*exp(EigValues_vec(r)*t_row)+A_r*v_3*exp(EigValues_vec(r)*t_row);
    
    %For Display only
    w_d_r=abs(imag(EigValues_vec(r)));
    zeta_r_w_r=-real(EigValues_vec(r));
    if w_d_r~=0 && mod(r,2)~=0  %Odd r
        disp('Free Response parameters');
        r %#ok<NOPRT>
         2*((-zeta_r_w_r*real(A_r)-w_d_r*imag(A_r))*v_1+real(A_r)*v_3)    %#ok<NOPRT> %cos coeff
        -2*((w_d_r*real(A_r)-zeta_r_w_r*imag(A_r))*v_1+imag(A_r)*v_3)     %#ok<NOPRT> %sin coeff
    end
end

x_rows(imag(x_rows)<100*eps)=real(x_rows(imag(x_rows)<100*eps));