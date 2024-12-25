function MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,...
                   isProportional,maxPhaseLag,display_EVD_Details)  %Optional arguments

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

[EigVectors_Normalized,s_q_vec]=MDOF_Eig_Visc(M_mat,C_mat,K_mat,isProportional,display_EVD_Details);

[w_p_vec,zeta_p_vec,w_d_p_vec]=pole2modal_visc(s_q_vec) %#ok<NOPRT>

%TF
H_s_mat=MDOF_TF_Visc(s_q_vec,EigVectors_Normalized);

%FRF & IRF labels
h_cols_Y_label_col=strings(n_RF_curves,1);
FRF_legend_str_vec=strings(n_RF_curves,1);
for ii=1:n_RF_curves
    h_cols_Y_label_col(ii)="$h_{"+m_row(ii)+','+n_row(ii)+'}(t)$';
    FRF_legend_str_vec(ii)="$H_{"+m_row(ii)+','+n_row(ii)+'}(f)$';
end

%FRF
H_w_n_m_cols=MDOF_FRF_Visc(s_q_vec,EigVectors_Normalized,w_column,m_row,n_row);
figure
ax_mag=plot_FRF_mag_phase(f_column,H_w_n_m_cols,false,[],[],[],[],[],maxPhaseLag);
legend(ax_mag,FRF_legend_str_vec,'interpreter','latex')
figure;
plot_FRF_Nyq(H_w_n_m_cols);
legend(FRF_legend_str_vec,'interpreter','latex')

%Antiresonance and minimum FRF
[w_11_AR,H_11_AR]=fminbnd(@(w) abs(MDOF_FRF_Visc(s_q_vec,EigVectors_Normalized,w,1,1)),w_d_p_vec(1),w_d_p_vec(2),optimset('TolX',1e-10))
[w_12_min,H_12_min]=fminbnd(@(w) abs(MDOF_FRF_Visc(s_q_vec,EigVectors_Normalized,w,2,1)),w_d_p_vec(1),w_d_p_vec(2),optimset('TolX',1e-10))
[w_22_AR,H_22_AR]=fminbnd(@(w) abs(MDOF_FRF_Visc(s_q_vec,EigVectors_Normalized,w,2,2)),w_d_p_vec(1),w_d_p_vec(2),optimset('TolX',1e-10))

%IRF
h_cols=MDOF_IRF_Visc(s_q_vec,EigVectors_Normalized,t_row.',m_row,n_row);
figure
tiledlayout(n_RF_curves,1,"TileSpacing","tight")
for n=1:n_RF_curves
    nexttile
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
x_str_col="x_{"+(1:N).'+'}(t)';

%Free response
x_rows=MDOF_Free_Response_Visc(M_mat,C_mat,s_q_vec,EigVectors_Normalized,x_0_col,x_dot_0_col,t_row);
figure
tiledlayout(N,1,"TileSpacing","tight")
for n=1:N
    nexttile
    plot(t_row,x_rows(n,:))
    ylabel("$"+x_str_col(n)+'$','interpreter','latex')
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
f_str="f_{1} (t)";
F_0_col=zeros(N,1);
F_0_col(1)=1;
w_F1=[0.5,0.9,1,1.1,w_11_AR/w_p_vec(1),w_12_min/w_p_vec(1)]*w_p_vec(1);
f_rows=zeros(N,n_points);
x_rows=zeros(N,n_points,2);
r_str_col=strings(N,1);
ignoreTransientVector=false;
if all(all(abs(C_mat)<=10000*eps))
    ignoreTransientVector=[true,ignoreTransientVector];
end
f_title_str="$"+f_str+"=\sin\left(\Omega_{1} t\right)$";
for ignoreTransient=ignoreTransientVector
    if ignoreTransient
        x_modified_str_col=strrep(x_str_col,'(','^{\mathrm{ss}}(');
        sameScale_y1_Vector=true;
    else
        x_modified_str_col=x_str_col;
        sameScale_y1_Vector=[false,true];
    end
    x_title_str_col="$"+x_modified_str_col+'$';
    figureTitle_col=["";""];
    if  all(abs(C_mat)<=10000*eps,"all") && ignoreTransient
        figureTitle_col="For undamped system, "+x_title_str_col+' never coincides with $'+x_str_col+'$';
    end

    for sameScale_y1=sameScale_y1_Vector
        for ii=1:length(w_F1)
            if w_F1(ii)==w_p_vec(1)
                r_str_col(ii)="$\Omega_{1}=\omega_{1}$";
            else
                r_str_col(ii)="$\Omega_{1}="+(w_F1(ii)/w_p_vec(1))+'\omega_{1}$';
            end
            w_F_col=zeros(N,1);
            w_F_col(1)=w_F1(ii);
            f_rows(ii,:)=F_0_col(1)*sin(w_F_col(1)*t_row);
            x_rows_temp=MDOF_Harmonic_Response_Visc(s_q_vec,EigVectors_Normalized,F_0_col,w_F_col,t_row,ignoreTransient);
            x_rows(ii,:,1)=x_rows_temp(1,:);
            x_rows(ii,:,2)=x_rows_temp(2,:);
        end
        r_str_col(end-1)="$\Omega_{1}=\omega_{1,1}^{\mathrm{AR}}$";
        r_str_col(end)="$\Omega_{1}=\omega_{1,2}^{\min}$";
        for nnn=1:2
            figure
            plot_Forced_Response_Vertically(t_row,x_rows(:,:,nnn),f_rows,figureTitle_col(nnn),f_title_str,r_str_col,x_title_str_col(nnn),sameScale_y1)
        end
    end
end

%Harmonic response 2
f_str="f_{2} (t)";
F_0_col=zeros(N,1);
F_0_col(2)=1;
w_F2=[w_12_min/w_p_vec(2),w_22_AR/w_p_vec(2),0.95,1,1.05,1.5]*w_p_vec(2);
f_title_str="$"+f_str+"=\sin\left(\Omega_{2} t\right)$";
for ignoreTransient=ignoreTransientVector
    if ignoreTransient
        x_modified_str_col=strrep(x_str_col,'(','^{\mathrm{ss}}(');
        sameScale_y1_Vector=true;
    else
        x_modified_str_col=x_str_col;
        sameScale_y1_Vector=[false,true];
    end
    x_title_str_col="$"+x_modified_str_col+'$';
    figureTitle_col=["";""];
    if  all(abs(C_mat)<=10000*eps,"all") && ignoreTransient
        figureTitle_col="For undamped system, "+x_title_str_col+' never coincides with $'+x_str_col+'$';
    end

    for sameScale_y1=sameScale_y1_Vector
        for ii=1:length(w_F2)
            if w_F2(ii)==w_p_vec(2)
                r_str_col(ii)="$\Omega_{2}=\omega_{2}$";
            else
                r_str_col(ii)="$\Omega_{2}="+(w_F2(ii)/w_p_vec(2))+'\omega_{2}$';
            end
            w_F_col=zeros(N,1);
            w_F_col(2)=w_F2(ii);
            f_rows(ii,:)=F_0_col(2)*sin(w_F_col(2)*t_row);
            x_rows_temp=MDOF_Harmonic_Response_Visc(s_q_vec,EigVectors_Normalized,F_0_col,w_F_col,t_row,ignoreTransient);
            x_rows(ii,:,1)=x_rows_temp(1,:);
            x_rows(ii,:,2)=x_rows_temp(2,:);
        end
        r_str_col(1)="$\Omega_{2}=\omega_{1,2}^{\min}$";
        r_str_col(2)="$\Omega_{2}=\omega_{2,2}^{\mathrm{AR}}$";
        for nnn=1:2
            figure
            plot_Forced_Response_Vertically(t_row,x_rows(:,:,nnn),f_rows,figureTitle_col(nnn),f_title_str,r_str_col,x_title_str_col(nnn),sameScale_y1)
        end
    end
end