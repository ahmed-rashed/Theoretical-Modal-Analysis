function [EigVectors_Normalized,s_vec]=MDOF_Eig_Visc(M_mat,C_mat,K_mat, ...
                 isPropotional,displayDetails)     %Optional arguments
if nargin<4
    isPropotional=false;
end

if nargin<5
    displayDetails=false;
end

P=size(M_mat,1);
if isPropotional || all(all(C_mat==0))    %Undamped or proportional
    [EigVectors_U,s_U_vec]=eig(K_mat,M_mat,'chol','vector');

    %Sort eigenvalues and corresponding eignvectors
    [~,Index]=sort(abs(s_U_vec));
    s_U_vec=s_U_vec(Index);
    EigVectors_U=EigVectors_U(:,Index);

    %Check the accuracy of eigendecomposition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s_U_prec=eps*P*max(abs(s_U_vec));

    if any(s_U_vec<=-s_U_prec)
        s_U_vec
        error('The matrix pencil (K,M) must be positive semi definite. That is; eigenvalues must be >=0')
    end

    if displayDetails
        s_U_vec

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Only necesary for better display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii=1:P
            EigVectors_U(:,ii)=EigVectors_U(:,ii)/EigVectors_U(1,ii);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EigVectors_U
    end

    IndexTemp=find(abs(s_U_vec)<=s_U_prec);
    if ~isempty(IndexTemp)
        if ~displayDetails
            s_U_vec
            EigVectors_U
        end
        warning('Calculated eigenvalues are inaccurate.')
        warning('Small eigenvalues and small elements in eigenvectors are manually reset to zero as follows:')
        disp('Press any key to continue');
        pause

        s_U_vec(IndexTemp)=0;s_U_vec
        EigVectors_U(abs(EigVectors_U)<=s_U_prec)=0
    end

    M_p_mat=EigVectors_U.'*M_mat*EigVectors_U;M_p_mat(abs(M_p_mat)<100*eps)=0;
    M_p_col=diag(M_p_mat);
    w_U_p_col=sqrt(s_U_vec);
    C_p_mat=EigVectors_U.'*C_mat*EigVectors_U;C_p_mat(abs(C_p_mat)<100*eps)=0;
    C_p_col=diag(C_p_mat);

    if displayDetails
        M_p_mat
        C_p_mat
    end

    %zeta_p_col=C_p_col/2./M_p_col./w_U_p_col;
    w_d_p_col=sqrt(w_U_p_col.^2-(C_p_col/2./M_p_col).^2);    %This is instead "w_U_p_mat.*sqrt(1-zeta_p_mat.^2)" to avoid the 0*inf in case w_U_p=0

    s_vec_temp1=-C_p_col/2./M_p_col-1i*w_d_p_col;
    s_vec_temp2=-C_p_col/2./M_p_col+1i*w_d_p_col;   %For overdamped proportional damping, Eigenvalues become real distinct

    s_vec=zeros(2*P,1);
    s_vec(1:2:2*P-1)=s_vec_temp1;
    s_vec(2:2:2*P)=s_vec_temp2;

    EigVectors_Normalized=zeros(P,2*P);
    EigVectors_Normalized(:,1:2:2*P-1)=EigVectors_U/sqrt(-1i*2*diag(w_d_p_col).*M_p_mat);
    EigVectors_Normalized(:,2:2:2*P)  =EigVectors_U/sqrt( 1i*2*diag(w_d_p_col).*M_p_mat);      %w_d_p_col may be complex for over damped modes
else    %Non-proportional
    [EigVectors_Normalized,s_vec]=quad_eig(K_mat,C_mat,M_mat);

    %Sort eigenvalues and corresponding eignvectors
    [~,Index]=sort(abs(imag(s_vec)));
    s_vec=s_vec(Index);
    EigVectors_Normalized=EigVectors_Normalized(:,Index);
end

%Only necesary for display
if displayDetails
    s_vec
    EigVectors_Normalized,

    for p=1:2*P
        if imag(s_vec(p))~=0 && mod(p,2)==0   %complex eigenvalue and even p
            continue
        end

        p
        if all(all(C_mat==0))    %Undamped
            A_p=EigVectors_U(:,(p+1)/2)*EigVectors_U(:,(p+1)/2).'/M_p_mat((p+1)/2,(p+1)/2)
        else
            A_p=EigVectors_Normalized(:,p)*EigVectors_Normalized(:,p).'
        end
    end
end
