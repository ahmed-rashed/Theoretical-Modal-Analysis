function H_s_mat=MDOF_TF_Visc(EigValues_vec,EigVectors_Normalized)

P=size(EigVectors_Normalized,1);

num=num2cell(zeros(P,P));
den=1;
H_s_mat=tf(num,den);
H_s_mat1=H_s_mat;
H_s_mat_SDOF=H_s_mat;
for q=1:2*P
    A_q=EigVectors_Normalized(:,q)*EigVectors_Normalized(:,q).';
    num=num2cell(A_q);
    den=[1,-EigValues_vec(q)];
    tf_i=tf(num,den);
    H_s_mat=H_s_mat+tf_i;
    H_s_mat_SDOF=H_s_mat_SDOF+tf_i;
    
    if mod(q,2)==0   %Even q
        %For Display only
        %%%%%%%%%%%%%%%%%%
        q
        H_s_mat_SDOF
        %%%%%%%%%%%%%%%%%%

        H_s_mat1=H_s_mat1+H_s_mat_SDOF;
        H_s_mat_SDOF=tf(num2cell(zeros(P,P)),1);
    end
end

