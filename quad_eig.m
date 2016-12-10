function [Epsi_normalized,Val_mat]=quad_eig(K,C,M)

L=2;    %Linearization type; 1 or 2

n=size(M,1);

if L==1 %L1 form
    NN=-K;
    A=[0*M,NN;-K,-C];
    B=[NN,0*M;0*M,M];
elseif L==2 %L2 form
    NN=M;
    A=[-K,0*M;0*M,NN];
    B=[C,M;NN,0*M];
end

[Phi,Val_mat] = eig(A,B);

%Check the accuracy of eigendecomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigValues_prec=eps*2*n*max(abs(diag(Val_mat)));

%     if any(diag(Val_mat)>=EigValues_prec)
%         Val_mat,
%         error('The matrix pencil (-K,M) must be negative semi definite. That is; eigenvalues must be >=0')
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