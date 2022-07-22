function MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,...
				   isProportional,maxPhaseLag,display_EVD_Details)	%Optional arguments

if nargin<10
    isProportional=false;
end

if nargin<11
    maxPhaseLag=[];
end

if nargin<12
    display_EVD_Details=false;
end

N=size(M_mat,1);
n_points=length(t_row);
n_f_points=length(f_column);
n_RF_curves=length(m_row);

w_column=2*pi*f_column;

[EigVectors_Normalized,EigValues_vec]=MDOF_Eig_Visc(M_mat,C_mat,K_mat,isProportional,display_EVD_Details);

[w_r_vec,zeta_r_vec,w_d_r_vec]=pole2modal_visc(EigValues_vec)

%TF
H_s_mat=MDOF_TF_Visc(EigValues_vec,EigVectors_Normalized);

%FRF & IRF labels
h_cols_Y_label_col=strings(n_RF_curves,1);
FRF_legend_str=strings(n_RF_curves,1);
for ii=1:n_RF_curves
    h_cols_Y_label_col(ii)="$h_{"+m_row(ii)+','+n_row(ii)+'}(t)$';
    FRF_legend_str(ii)="$H_{"+m_row(ii)+','+n_row(ii)+'}(f)$';
end

%FRF
H_w_n_m_cols=MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,w_column,m_row,n_row);
figure
ax_mag=plot_FRF_mag_phase(f_column,H_w_n_m_cols,false,[],[],[],[],[],maxPhaseLag);
legend(ax_mag,FRF_legend_str,'interpreter','latex')
figure;
plot_FRF_Nyq(H_w_n_m_cols);
legend(FRF_legend_str,'interpreter','latex')

%Antiresonance and minimum FRF
[w_11_AR,H_11_AR]=fminbnd(@(w) abs(MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,w,1,1)),w_d_r_vec(1),w_d_r_vec(2),optimset('TolX',1e-10))
[w_12_min,H_12_min]=fminbnd(@(w) abs(MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,w,2,1)),w_d_r_vec(1),w_d_r_vec(2),optimset('TolX',1e-10))
[w_22_AR,H_22_AR]=fminbnd(@(w) abs(MDOF_FRF_Visc(EigValues_vec,EigVectors_Normalized,w,2,2)),w_d_r_vec(1),w_d_r_vec(2),optimset('TolX',1e-10))

%IRF
h_cols=MDOF_IRF_Visc(EigValues_vec,EigVectors_Normalized,t_row.',m_row,n_row);
figure
for n=1:n_RF_curves
    subplot(n_RF_curves,1,n)
    plot(t_row,h_cols(:,n).')
    ylabel(h_cols_Y_label_col(n),'interpreter','latex')
    if n==1
        title('IRF','interpreter','latex')
    end
    if n==n_RF_curves
        xlabel('$t$','interpreter','latex')
    else
        set(gca,'XTickLabel',[]);
    end
end

%Response Labels
x_ylabel_col=strings(N,1);
for ii=1:N
    x_ylabel_col(ii)="$x_{"+ii+'}(t)$';
end
x_new_ylabel_col=x_ylabel_col;

%Free response
x_rows=MDOF_Free_Response_Visc(M_mat,C_mat,EigValues_vec,EigVectors_Normalized,x_0_col,x_dot_0_col,t_row);
figure
for n=1:N
    subplot(N,1,n)
    plot(t_row,x_rows(n,:))
    ylabel(x_ylabel_col(n),'interpreter','latex')
    if n==1
        title('Free response','interpreter','latex')
    end
    if n==N
        xlabel('$t$','interpreter','latex')
    else
        set(gca,'XTickLabel',[]);
    end
end

%Harmonic response 1
F_0_col=zeros(N,1);
F_0_col(1)=1;
w_F1=[0.5,0.9,1,1.1,w_11_AR/w_r_vec(1),w_12_min/w_r_vec(1)]*w_r_vec(1);
f_rows=zeros(N,n_points);
x_rows1=zeros(N,n_points);
x_rows2=zeros(N,n_points);
f_rows_labels_col=strings(N,1);
ignoreTransientVector=false;
if all(all(abs(C_mat)<=10000*eps))
    ignoreTransientVector=[true,ignoreTransientVector];
end
for ignoreTransient=ignoreTransientVector
    if ignoreTransient
        x_new_ylabel_col(1)=strrep(x_ylabel_col(1),'(','^{\mathrm{ss}}(');
        x_new_ylabel_col(2)=strrep(x_ylabel_col(2),'(','^{\mathrm{ss}}(');
        sameScale_y1_Vector=true;
    else
        x_new_ylabel_col=x_ylabel_col;
        sameScale_y1_Vector=[false,true];
    end
    figureTitle1=x_new_ylabel_col(1)+' due to $f_{1} (t)=\sin\left(\Omega_{1}t\right)$';
    figureTitle2=x_new_ylabel_col(2)+' due to $f_{1} (t)=\sin\left(\Omega_{1}t\right)$';
    if  all(all(abs(C_mat)<=10000*eps))
        figureTitle1=figureTitle1+' for undamped system';
        figureTitle2=figureTitle2+' for undamped system';
    end
        
    for sameScale_y1=sameScale_y1_Vector
        if  ignoreTransient && all(all(abs(C_mat)<=10000*eps))
            figureTitle1=figureTitle1+' \underline{(never coincides with '+x_ylabel_col(1)+', but matches $H_{1,1}(\omega)$)}';
            figureTitle2=figureTitle2+' \underline{(never coincides with '+x_ylabel_col(2)+', but matches $H_{2,1}(\omega)$)}';
        end

        for ii=1:length(w_F1)
            if w_F1(ii)==w_r_vec(1)
                f_rows_labels_col(ii)="$f_{1} (t),:\Omega_{1}=\omega_{1}$";
            else
                f_rows_labels_col(ii)="$f_{1} (t),:\Omega_{1}="+(w_F1(ii)/w_r_vec(1))+'\omega_{1}$';
            end
            w_F_col=zeros(N,1);
            w_F_col(1)=w_F1(ii);
            f_rows(ii,:)=F_0_col(1)*sin(w_F_col(1)*t_row);
            x_rows_temp=MDOF_Harmonic_Response_Visc(EigValues_vec,EigVectors_Normalized,F_0_col,w_F_col,t_row,ignoreTransient);
            x_rows1(ii,:)=x_rows_temp(1,:);
            x_rows2(ii,:)=x_rows_temp(2,:);
        end
        f_rows_labels_col(end-1)="$f_{1} (t),:\Omega_{1}=\omega_{1,1}^{\mathrm{AR}}$";
        f_rows_labels_col(end)="$f_{1} (t),:\Omega_{1}=\omega_{1,2}^{\min}$";
        figure
        plot_Forced_Response_Vertically(t_row,x_rows1,x_new_ylabel_col(1),f_rows,f_rows_labels_col,figureTitle1,sameScale_y1)
        figure
        plot_Forced_Response_Vertically(t_row,x_rows2,x_new_ylabel_col(2),f_rows,f_rows_labels_col,figureTitle2,sameScale_y1)
    end
end

%Harmonic response 2
F_0_col=zeros(N,1);
F_0_col(2)=1;
w_F2=[w_12_min/w_r_vec(2),w_22_AR/w_r_vec(2),0.95,1,1.05,1.5]*w_r_vec(2);
for ignoreTransient=ignoreTransientVector
    if ignoreTransient
        x_new_ylabel_col(1)=strrep(x_ylabel_col(1),'(','^{\mathrm{ss}}(');
        x_new_ylabel_col(2)=strrep(x_ylabel_col(2),'(','^{\mathrm{ss}}(');
        sameScale_y1_Vector=true;
    else
        x_new_ylabel_col=x_ylabel_col;
        sameScale_y1_Vector=[false,true];
    end
    figureTitle1=x_new_ylabel_col(1)+' due to $f_{2} (t)=\sin\left(\Omega_{2}t\right)$';
    figureTitle2=x_new_ylabel_col(2)+' due to $f_{2} (t)=\sin\left(\Omega_{2}t\right)$';
    if  all(all(abs(C_mat)<=10000*eps))
        figureTitle1=figureTitle1+' for undamped system';
        figureTitle2=figureTitle2+' for undamped system';
    end
    
    for sameScale_y1=sameScale_y1_Vector
        if  ignoreTransient && all(all(abs(C_mat)<=10000*eps))
            figureTitle1=figureTitle1+' \underline{(never coincides with '+x_ylabel_col(1)+', but matches $H_{1,2}(\omega)$)}';
            figureTitle2=figureTitle2+' \underline{(never coincides with '+x_ylabel_col(2)+', but matches $H_{2,2}(\omega)$)}';
        end

        for ii=1:length(w_F2)
            f_rows_labels_col(ii)="$f_{2} (t),:\Omega_{2}="+(w_F2(ii)/w_r_vec(2))+'\omega_{2}$';
            w_F_col=zeros(N,1);
            w_F_col(2)=w_F2(ii);
            f_rows(ii,:)=F_0_col(2)*sin(w_F_col(2)*t_row);
            x_rows_temp=MDOF_Harmonic_Response_Visc(EigValues_vec,EigVectors_Normalized,F_0_col,w_F_col,t_row,ignoreTransient);
            x_rows1(ii,:)=x_rows_temp(1,:);
            x_rows2(ii,:)=x_rows_temp(2,:);
        end
        f_rows_labels_col(1)="$f_{2} (t),:\Omega_{2}=\omega_{1,2}^{\min}$";
        f_rows_labels_col(2)="$f_{2} (t),:\Omega_{2}=\omega_{2,2}^{\mathrm{AR}}$";
        figure
        plot_Forced_Response_Vertically(t_row,x_rows1,x_new_ylabel_col(1),f_rows,f_rows_labels_col,figureTitle1,sameScale_y1)
        figure
        plot_Forced_Response_Vertically(t_row,x_rows2,x_new_ylabel_col(2),f_rows,f_rows_labels_col,figureTitle2,sameScale_y1)
    end
end
