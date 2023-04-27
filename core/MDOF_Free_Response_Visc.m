function x_rows= ...
MDOF_Free_Response_Visc(M,C,s_q_vec,EigVectors_Normalized,x_0_col,x_dot_0_col,t_row)

P=size(EigVectors_Normalized,1);
n_col=size(t_row,2);
if size(x_0_col,1)~=P,error('initial conditions must have the same length as the order of the square matrix M '),end

x_rows=zeros(P,n_col);
v_1=M*x_0_col;
v_3=C*x_0_col+M*x_dot_0_col;
for p=1:2*P
    A_p=EigVectors_Normalized(:,p)*EigVectors_Normalized(:,p).';
    x_rows=x_rows+A_p*(s_q_vec(p)*v_1*exp(s_q_vec(p)*t_row)+v_3*exp(s_q_vec(p)*t_row));
    
    %For Display only
    w_d_p=abs(imag(s_q_vec(p)));
    zeta_p_w_p=-real(s_q_vec(p));
    if w_d_p~=0 && mod(p,2)~=0  %Odd p
        disp('Free Response parameters');
        p %#ok<NOPRT>
         2*((-zeta_p_w_p*real(A_p)-w_d_p*imag(A_p))*v_1+real(A_p)*v_3)    %#ok<NOPRT> %cos coeff
        -2*((w_d_p*real(A_p)-zeta_p_w_p*imag(A_p))*v_1+imag(A_p)*v_3)     %#ok<NOPRT> %sin coeff
    end
end

mask=imag(x_rows)<100*eps;
x_rows(mask)=real(x_rows(mask));