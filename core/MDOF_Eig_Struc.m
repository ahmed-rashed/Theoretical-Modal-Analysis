function [EigVectors_Normalized, EigValues_vec]= ...
MDOF_Eig_Struc(M_mat, D_mat, K_mat, ...
				 displayDetails)     %Optional arguments

if nargin<4
    displayDetails=false;
end

N=size(M_mat,1);

[EigVectors_H,EigValues_H_mat]=eig(-(K_mat+1i*D_mat),M_mat);
EigValues_H_vec=diag(EigValues_H_mat);

%Sort eigenvalues and corresponding eignvectors
[~,Index]=sort(real(EigValues_H_vec));
EigValues_H_vec=EigValues_H_vec(Index);
EigVectors_H=EigVectors_H(:,Index);

%Check the accuracy of eigendecomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigValues_prec=eps*N*max(abs(EigValues_H_vec));

if any(real(EigValues_H_vec)>=EigValues_prec)
    EigValues_H_vec
    error('The matrix pencil (-(K+1i*D),M) eigenvalues must have real part <=0')
end

EigValues_H_vec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Only for better display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:N
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

M_r_mat=EigVectors_H.'*M_mat*EigVectors_H;M_r_mat(abs(M_r_mat)<100*eps)=0;
M_r_col=diag(M_r_mat);
w_U_r_col=sqrt(-EigValues_H_vec);
C_r_mat=EigVectors_H.'*D_mat*EigVectors_H;C_r_mat(abs(C_r_mat)<100*eps)=0;
C_r_col=diag(C_r_mat);

if displayDetails
    M_r_mat
    C_r_mat
end

%zeta_r_col=C_r_col/2./M_r_col./w_U_r_col;
w_d_r_col=sqrt(w_U_r_col.^2-(C_r_col/2./M_r_col).^2);    %This is instead "w_U_r_mat.*sqrt(1-zeta_r_mat.^2)" to avoid the 0*inf in case w_U_r=0

EigValues_vec_temp1=-C_r_col/2./M_r_col-1i*w_d_r_col;
EigValues_vec_temp2=-C_r_col/2./M_r_col+1i*w_d_r_col;   %Eigenvalues not necessarily complex conjugate pairs

EigValues_vec=zeros(2*N,1);
EigValues_vec(1:2:2*N-1)=EigValues_vec_temp1;
EigValues_vec(2:2:2*N)=EigValues_vec_temp2;

EigVectors_Normalized=zeros(N,2*N);
EigVectors_Normalized(:,1:2:2*N-1)=EigVectors_H/sqrt(-i*2*diag(w_d_r_col).*M_r_mat);
EigVectors_Normalized(:,2:2:2*N)=EigVectors_H/sqrt(i*2*diag(w_d_r_col).*M_r_mat);      %w_d_r_col may be complex for over damped modes

%Only necesary for display
if displayDetails
    EigValues_vec
    EigVectors_Normalized,
    
    for r=1:2*N
        if imag(EigValues_vec(r))~=0 && mod(r,2)==0   %complex eigenvalue and even r
            continue
        end

        r
        if all(all(D_mat==0))    %Undamped
            A_r=EigVectors_H(:,(r+1)/2)*EigVectors_H(:,(r+1)/2).'/M_r_mat((r+1)/2,(r+1)/2)
        else
            A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).'
        end
    end
end
