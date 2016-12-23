function [EigVectors_Normalized, EigValues_mat]= ...
MDOF_Eig_Visc(M, C, K, ...
				 isPropotional)     %Optional arguments
if nargin<4
    isPropotional=false;
end

N=size(M,1);
if isPropotional || all(all(C==0))    %Undamped or proportional
    [EigVectors_U,EigValues_U_mat]=eig(-K,M);
    
%Only necesary for display
for ii=1:N
    EigVectors_U(:,ii)=EigVectors_U(:,ii)/EigVectors_U(1,ii);
end
[Temp_Val,Index]=sort(abs(diag(EigValues_U_mat)));
EigValues_U_mat=EigValues_U_mat(Index,Index);
EigVectors_U=EigVectors_U(:,Index);
EigValues_U_mat,EigVectors_U
    
    M_r=EigVectors_U.'*M*EigVectors_U
    w_U_r=sqrt(-EigValues_U_mat)

    C_r=EigVectors_U.'*C*EigVectors_U
    zeta_r=C_r/(2*M_r*w_U_r)
    w_d_r=w_U_r.*sqrt(1-zeta_r.^2)
    EigValues_vec_temp=diag(-zeta_r*w_U_r+i*w_d_r);

    EigValues_vec=zeros(2*N,1);
    EigValues_vec(1:2:2*N-1)=EigValues_vec_temp;
    EigValues_vec(2:2:2*N)=conj(EigValues_vec_temp);
    EigValues_mat=diag(EigValues_vec)
    
    EigVectors_Normalized=zeros(N,2*N);
    EigVectors_Normalized(:,1:2:2*N-1)=EigVectors_U/sqrt(i*2*w_d_r.*M_r);
    EigVectors_Normalized(:,2:2:2*N)=conj(EigVectors_Normalized(:,1:2:2*N-1))
else    %Non-proportional
    [EigVectors_Normalized,EigValues_mat]=quad_eig(K,C,M);
end

[Temp_Val,Index]=sort(abs(imag(diag(EigValues_mat))));
EigValues_mat=EigValues_mat(Index,Index);
EigVectors_Normalized=EigVectors_Normalized(:,Index);

% %Only necesary for display
% for r=1:2:2*N
%     r
%     if all(all(C==0))    %Undamped
%         A_r=EigVectors_U(:,(r+1)/2)*EigVectors_U(:,(r+1)/2).'/M_r((r+1)/2,(r+1)/2)
%     else
%         A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).'
%     end
% end
