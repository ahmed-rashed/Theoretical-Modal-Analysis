function [EigVectors_Normalized, EigValues_vec]= ...
MDOF_Eig_Visc(M_mat, C_mat, K_mat, ...
				 isPropotional,displayDetails)     %Optional arguments
if nargin<4
    isPropotional=false;
end

if nargin<5
    displayDetails=false;
end

N=size(M_mat,1);
if isPropotional || all(all(C_mat==0))    %Undamped or proportional
    [EigVectors_U,EigValues_U_mat]=eig(-K_mat,M_mat);
    EigValues_U_vec=diag(EigValues_U_mat);
    
    %Sort eigenvalues and corresponding eignvectors
    [~,Index]=sort(abs(EigValues_U_vec));
    EigValues_U_vec=EigValues_U_vec(Index);
    EigVectors_U=EigVectors_U(:,Index);

    %Check the accuracy of eigendecomposition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EigValues_prec=eps*N*max(abs(EigValues_U_vec));

    if any(EigValues_U_vec>=EigValues_prec)
        EigValues_U_vec
        error('The matrix pencil (-K,M) must be negative semi definite. That is; eigenvalues must be >=0')
    end

    EigValues_U_vec
    EigVectors_U_temp=EigVectors_U;
    for ii=1:N
        EigVectors_U_temp(:,ii)=EigVectors_U_temp(:,ii)/EigVectors_U_temp(1,ii);
    end
    EigVectors_U_temp
    IndexTemp=find(abs(EigValues_U_vec)<=EigValues_prec);
    if ~isempty(IndexTemp)
        warning('Calculated eigenvalues are inaccurate.')
        warning('Small eigenvalues and small elements in eigenvectors are manually reset to zero as follows:')
        disp('Press any key to continue');
        pause
        
        EigValues_U_vec(IndexTemp)=0;EigValues_U_vec
        EigVectors_U(abs(EigVectors_U)<EigValues_prec)=0
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Only necesary for better display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:N
        EigVectors_U(:,ii)=EigVectors_U(:,ii)/EigVectors_U(1,ii);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    M_r_mat=EigVectors_U.'*M_mat*EigVectors_U;M_r_mat(abs(M_r_mat)<100*eps)=0;
    M_r_col=diag(M_r_mat);
    w_U_r_col=sqrt(-EigValues_U_vec);
    C_r_mat=EigVectors_U.'*C_mat*EigVectors_U;C_r_mat(abs(C_r_mat)<100*eps)=0;
    C_r_col=diag(C_r_mat);
    
    if displayDetails
        M_r_mat
        C_r_mat
    end
    
    %zeta_r_col=C_r_col/2./M_r_col./w_U_r_col;
    w_d_r_col=sqrt(w_U_r_col.^2-(C_r_col/2./M_r_col).^2);    %This is instead "w_U_r_mat.*sqrt(1-zeta_r_mat.^2)" to avoid the 0*inf in case w_U_r = 0
    
    EigValues_vec_temp1=-C_r_col/2./M_r_col-1i*w_d_r_col;
    EigValues_vec_temp2=-C_r_col/2./M_r_col+1i*w_d_r_col;   %Eigenvalues not necessarily complex conjugate pairs

    EigValues_vec=zeros(2*N,1);
    EigValues_vec(1:2:2*N-1)=EigValues_vec_temp1;
    EigValues_vec(2:2:2*N)=EigValues_vec_temp2;
    
    EigVectors_Normalized=zeros(N,2*N);
    EigVectors_Normalized(:,1:2:2*N-1)=EigVectors_U/sqrt(-i*2*diag(w_d_r_col).*M_r_mat);
    EigVectors_Normalized(:,2:2:2*N)=EigVectors_U/sqrt(i*2*diag(w_d_r_col).*M_r_mat);      %w_d_r_col may be complex for over damped modes
else    %Non-proportional
    [EigVectors_Normalized,EigValues_vec]=quad_eig(K_mat,C_mat,M_mat);
    
    %Sort eigenvalues and corresponding eignvectors
    [~,Index]=sort(abs(imag(EigValues_vec)));
    EigValues_vec=EigValues_vec(Index);
    EigVectors_Normalized=EigVectors_Normalized(:,Index);
end

%Only necesary for display
if displayDetails
    EigValues_vec
    EigVectors_Normalized,
    
    for r=1:2*N
        if imag(EigValues_vec(r))~=0 && mod(r,2)==0   %complex eigenvalue and even r
            continue
        end

        r
        if all(all(C_mat==0))    %Undamped
            A_r=EigVectors_U(:,(r+1)/2)*EigVectors_U(:,(r+1)/2).'/M_r_mat((r+1)/2,(r+1)/2)
        else
            A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).'
        end
    end
end