function [Epsi_normalized,Val_vec]=quad_eig(K_mat,C_mat,M_mat)

L=2;    %Linearization type; 1 or 2

N=size(M_mat,1);

if L==1 %L1 form
    N_mat=-K_mat;
    A=[0*M_mat,N_mat;-K_mat,-C_mat];
    B=[N_mat,0*M_mat;0*M_mat,M_mat];
elseif L==2 %L2 form
    N_mat=M_mat;
    A=[-K_mat,0*M_mat;0*M_mat,N_mat];
    B=[C_mat,M_mat;N_mat,0*M_mat];
end

[Phi,Val_vec]=eig(A,B,"vector");

%Check the accuracy of eigendecomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigValues_prec=eps*2*N*max(abs(Val_vec));

%     if any(Val_vec>=EigValues_prec)
%         Val_vec,
%         error('The matrix pencil (-K_mat,M_mat) must be negative semi definite. That is; eigenvalues must be >=0')
%     end

IndexTemp=find(abs(Val_vec)<=EigValues_prec);
if ~isempty(IndexTemp)
    Val_vec,
    Phi,

    warning('Calculated eigenvalues are inaccurate.')
    warning('Small eigenvalues and small elements in eigenvectors are manually reset to zero as follows:')
    disp('Press any key to continue');pause

    Val_vec(IndexTemp)=0
    Phi(abs(Phi)<EigValues_prec)=0
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_r_vec=diag(Phi.'*B*Phi);

Phi_normalized=Phi;
if L==1 %L1
	for ii=1:2*N
		Phi_normalized(:,ii)=Phi_normalized(:,ii)*sqrt(Val_vec(ii)/B_r_vec(ii));
	end
elseif L==2 %L2
	for ii=1:2*N
		Phi_normalized(:,ii)=Phi_normalized(:,ii)/sqrt(B_r_vec(ii));
	end
end

Epsi_normalized=Phi_normalized(1:N,:);
