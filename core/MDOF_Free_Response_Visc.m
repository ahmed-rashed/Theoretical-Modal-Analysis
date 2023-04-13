function x_rows= ...
MDOF_Free_Response_Visc(M,C,EigValues_vec,EigVectors_Normalized,x_0_col,x_dot_0_col,t_row)

P=size(EigVectors_Normalized,1);
n_col=size(t_row,2);
if size(x_0_col,1)~=P,error('initial conditions must have the same length as the order of the square matrix M '),end

x_rows=zeros(P,n_col);
v_1=M*x_0_col;
v_3=C*x_0_col+M*x_dot_0_col;
for p=1:2*P
    A_r=EigVectors_Normalized(:,p)*EigVectors_Normalized(:,p).';
    x_rows=x_rows+A_r*(EigValues_vec(p)*v_1*exp(EigValues_vec(p)*t_row)+v_3*exp(EigValues_vec(p)*t_row));
    
    %For Display only
    w_d_r=abs(imag(EigValues_vec(p)));
    zeta_r_w_r=-real(EigValues_vec(p));
    if w_d_r~=0 && mod(p,2)~=0  %Odd p
        disp('Free Response parameters');
        p %#ok<NOPRT>
         2*((-zeta_r_w_r*real(A_r)-w_d_r*imag(A_r))*v_1+real(A_r)*v_3)    %#ok<NOPRT> %cos coeff
        -2*((w_d_r*real(A_r)-zeta_r_w_r*imag(A_r))*v_1+imag(A_r)*v_3)     %#ok<NOPRT> %sin coeff
    end
end

mask=imag(x_rows)<100*eps;
x_rows(mask)=real(x_rows(mask));