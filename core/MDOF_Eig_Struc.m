function [EigVectors_Normalized,EigValues_vec]= ...
MDOF_Eig_Struc(M_mat,D_mat,K_mat,...
				 displayDetails)     %Optional arguments

if nargin<4
    displayDetails=false;
end

P=size(M_mat,1);

[EigVectors_H,EigValues_H_mat]=eig(-(K_mat+1i*D_mat),M_mat);
EigValues_H_vec=diag(EigValues_H_mat);

%Sort eigenvalues and corresponding eignvectors
[~,Index]=sort(real(EigValues_H_vec));
EigValues_H_vec=EigValues_H_vec(Index);
EigVectors_H=EigVectors_H(:,Index);

%Check the accuracy of eigendecomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigValues_prec=eps*P*max(abs(EigValues_H_vec));

if any(real(EigValues_H_vec)>=EigValues_prec)
    EigValues_H_vec
    error('The matrix pencil (-(K+1i*D),M) eigenvalues must have real part <=0')
end

EigValues_H_vec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Only for better display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:P
    EigVectors_H(:,ii)=EigVectors_H(:,ii)/EigVectors_H(1,ii);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigVectors_H
IndexTemp=find(abs(EigValues_H_vec)<=EigValues_prec);
if ~isempty(IndexTemp)
    warning('Calculated eigenvalues are inaccurate.')
    warning('Small eigenvalues and small elements in eigenvectors are manually reset to zero as follows:')
    disp('Press any key to continue');
    pause

    EigValues_H_vec(IndexTemp)=0;EigValues_H_vec
    EigVectors_H(abs(EigVectors_H)<=EigValues_prec)=0
end

M_p_mat=EigVectors_H.'*M_mat*EigVectors_H;M_p_mat(abs(M_p_mat)<100*eps)=0;
M_p_col=diag(M_p_mat);
w_U_p_col=sqrt(-EigValues_H_vec);
C_p_mat=EigVectors_H.'*D_mat*EigVectors_H;C_p_mat(abs(C_p_mat)<100*eps)=0;
C_p_col=diag(C_p_mat);

if displayDetails
    M_p_mat
    C_p_mat
end

%zeta_p_col=C_p_col/2./M_p_col./w_U_p_col;
w_d_p_col=sqrt(w_U_p_col.^2-(C_p_col/2./M_p_col).^2);    %This is instead "w_U_p_mat.*sqrt(1-zeta_p_mat.^2)" to avoid the 0*inf in case w_U_p=0

EigValues_vec_temp1=-C_p_col/2./M_p_col-1i*w_d_p_col;
EigValues_vec_temp2=-C_p_col/2./M_p_col+1i*w_d_p_col;   %Eigenvalues not necessarily complex conjugate pairs

EigValues_vec=zeros(2*P,1);
EigValues_vec(1:2:2*P-1)=EigValues_vec_temp1;
EigValues_vec(2:2:2*P)=EigValues_vec_temp2;

EigVectors_Normalized=zeros(P,2*P);
EigVectors_Normalized(:,1:2:2*P-1)=EigVectors_H/sqrt(-i*2*diag(w_d_p_col).*M_p_mat);
EigVectors_Normalized(:,2:2:2*P)=EigVectors_H/sqrt(i*2*diag(w_d_p_col).*M_p_mat);      %w_d_p_col may be complex for over damped modes

%Only necesary for display
if displayDetails
    EigValues_vec
    EigVectors_Normalized,
    
    for p=1:2*P
        if imag(EigValues_vec(p))~=0 && mod(p,2)==0   %complex eigenvalue and even p
            continue
        end

        p
        if all(all(D_mat==0))    %Undamped
            A_p=EigVectors_H(:,(p+1)/2)*EigVectors_H(:,(p+1)/2).'/M_p_mat((p+1)/2,(p+1)/2)
        else
            A_p=EigVectors_Normalized(:,p)*EigVectors_Normalized(:,p).'
        end
    end
end
