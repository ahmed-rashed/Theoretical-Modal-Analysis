function H_s_mat=MDOF_TF_Visc(EigValues_vec,EigVectors_Normalized)

N=size(EigVectors_Normalized,1);

num=num2cell(zeros(N,N));
den=1;
H_s_mat=tf(num,den);
H_s_mat1=H_s_mat;
H_s_mat_SDOF=H_s_mat;
for r=1:2*N
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    num=num2cell(A_r);
    den=[1,-EigValues_vec(r)];
    tf_i=tf(num,den);
    H_s_mat=H_s_mat+tf_i;
    H_s_mat_SDOF=H_s_mat_SDOF+tf_i;
    
    if mod(r,2)==0   %Even r
        %For Display only
        %%%%%%%%%%%%%%%%%%
        r
        H_s_mat_SDOF
        %%%%%%%%%%%%%%%%%%

        H_s_mat1=H_s_mat1+H_s_mat_SDOF;
        H_s_mat_SDOF=tf(num2cell(zeros(N,N)),1);
    end
end

