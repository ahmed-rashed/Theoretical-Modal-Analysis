function H_s_mat=MDOF_TF_Visc(EigValues_vec, EigVectors_Normalized)

temp=abs(EigValues_vec);
if max(abs(temp(1:2:end)-temp(2:2:end)))>100*eps
    error('This function is designed for complex conjugate poles. Please upgrade this function for a general case poles');
end

N=size(EigVectors_Normalized,1);

H_s_mat=tf(zeros(N,N,1,1,2*N));
num=num2cell(zeros(N,N));
den=1;
H_s_mat_SDOF=tf(num,den);
for r=1:2*N
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    num=num2cell(A_r);
    den=[1,-EigValues_vec(r)];
    H_s_mat_SDOF=H_s_mat_SDOF+tf(num,den);
    
    if imag(EigValues_vec(r))~=0 && mod(r,2)~=0   %complex eigenvalue and odd r
        continue
    end

%For Display only
r
H_s_mat_SDOF

    H_s_mat(:,:,1,1,r)=H_s_mat_SDOF;
    H_s_mat_SDOF=tf(num2cell(zeros(N,N)),1);
end

