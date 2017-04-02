function [Epsi_normalized,Val_mat]=quad_eig(K_mat,C_mat,M_mat)

L=2;    %Linearization type; 1 or 2

n=size(M_mat,1);

if L==1 %L1 form
    N_mat=-K_mat;
    A=[0*M_mat,N_mat;-K_mat,-C_mat];
    B=[N_mat,0*M_mat;0*M_mat,M_mat];
elseif L==2 %L2 form
    N_mat=M_mat;
    A=[-K_mat,0*M_mat;0*M_mat,N_mat];
    B=[C_mat,M_mat;N_mat,0*M_mat];
end

[Phi,Val_mat] = eig(A,B);

%Check the accuracy of eigendecomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigValues_prec=eps*2*n*max(abs(diag(Val_mat)));

%     if any(diag(Val_mat)>=EigValues_prec)
%         Val_mat,
%         error('The matrix pencil (-K_mat,M_mat) must be negative semi definite. That is; eigenvalues must be >=0')
%     end

IndexTemp=find(abs(diag(Val_mat))<=EigValues_prec);
if ~isempty(IndexTemp)
    diag(Val_mat),
    Phi,

    warning('Calculated eigenvalues are inaccurate.')
    warning('Small eigenvalues and small elements in eigenvectors are manually reset to zero as follows:')
    disp('Press any key to continue');pause

    Val_mat(IndexTemp,IndexTemp)=0;diag(Val_mat)
    Phi(abs(Phi)<EigValues_prec)=0
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_r_mat=Phi.'*B*Phi;
B_r_mat(abs(B_r_mat)<100*eps)=0;

if L==1 %L1
    Phi_normalized=Phi*sqrt(Val_mat/B_r_mat);
elseif L==2 %L2
    Phi_normalized=Phi/sqrt(B_r_mat);
end

Epsi_normalized=Phi_normalized(1:n,:);