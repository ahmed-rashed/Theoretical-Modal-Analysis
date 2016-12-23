function MDOF_Prob(M,C,K,x_0_col,x_dot_0_col,n_row,m_row,t_row,f_column, ...
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

N=size(M,1);
n_points=length(t_row);
n_f_points=length(f_column);

w_column=2*pi*f_column;

[EigVectors_Normalized, EigValues_mat]=MDOF_Eig_Visc(M, C, K,isProportional,display_EVD_Details);

[w_r_vec, zeta_r_vec, w_d_r_vec]=MDOF_Modal_Param_Visc(EigValues_mat)

%TF
H_s_mat=MDOF_TF_Visc(EigValues_mat, EigVectors_Normalized);

%FRF & IRF labels
h_cols_Y_label_col=cell(length(n_row),1);
FRF_legend_str=cell(length(n_row),1);
for ii=1:length(n_row)
    h_cols_Y_label_col(ii)={['$h_{',int2str(n_row(ii)),',',int2str(m_row(ii)),'}$']};
    FRF_legend_str(ii)={['$H_{',int2str(n_row(ii)),',',int2str(m_row(ii)),'}$']};
end

%FRF
H_w_n_m_cols=MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, w_column, n_row, m_row);
figure
ax_mag=plot_FRF_mag_phase(f_column,H_w_n_m_cols,false,[],[],[],[],[],maxPhaseLag);
legend(ax_mag,FRF_legend_str,'interpreter','latex')
figure;
plot_FRF_Nyq(H_w_n_m_cols);
legend(FRF_legend_str,'interpreter','latex')

%Antiresonance and minimum FRF
[w_AR_11,H_11_AR]=fminbnd(@(w) abs(MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, w, 1, 1,false)), w_d_r_vec(1),w_d_r_vec(2),optimset('TolX',1e-10))
[w_min_12,H_12_min]=fminbnd(@(w) abs(MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, w, 2, 1,false)), w_d_r_vec(1),w_d_r_vec(2),optimset('TolX',1e-10))
[w_AR_22,H_22_AR]=fminbnd(@(w) abs(MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, w, 2, 2,false)), w_d_r_vec(1),w_d_r_vec(2),optimset('TolX',1e-10))

%IRF
h_cols=MDOF_IRF_Visc(EigValues_mat, EigVectors_Normalized, t_row.', n_row, m_row);
figure
plotResponse_ForceVertically(t_row,h_cols.',h_cols_Y_label_col,[],[],['IRF'])

%Response Labels
x_rows_Y_label_col=cell(N,1);
for ii=1:N
    x_rows_Y_label_col{ii}=['$x_{',int2str(ii),'}(t)$'];
end

%Free response
x_rows=MDOF_Free_Response_Visc(M,C, EigValues_mat, EigVectors_Normalized, x_0_col, x_dot_0_col, t_row);
figure
plotResponse_ForceVertically(t_row,x_rows,x_rows_Y_label_col,[],[],['Free response']);

%Harmonic response 1
F_0_col=zeros(N,1);
F_0_col(1)=1;
w_AR=sqrt(w_d_r_vec(1)^2+w_d_r_vec(2)^2)/sqrt(2);
w_F1=[0.5*w_d_r_vec(1),0.9*w_d_r_vec(1),w_d_r_vec(1),1.1*w_d_r_vec(1),w_AR_11,w_min_12];
f_rows=zeros(N,n_points);
x_rows1=zeros(N,n_points);
x_rows2=zeros(N,n_points);
f_rows_labels_col=cell(N,1);
ignoreTransientVector=false;
if all(all(abs(C)<=10000*eps))
    ignoreTransientVector=[true,ignoreTransientVector];
end
for ignoreTransient=ignoreTransientVector
    figureTitle1=[x_rows_Y_label_col{1},' due to $f_{1} (t)=\sin\left(w_{1}t\right)$'];
    figureTitle2=[x_rows_Y_label_col{2},' due to $f_{1} (t)=\sin\left(w_{1}t\right)$'];

    if ignoreTransient
        sameScale_y1_Vector=true;
    else
        sameScale_y1_Vector=[false,true];
    end
    
    for sameScale_y1=sameScale_y1_Vector
        if ignoreTransient
            figureTitle1=[figureTitle1,'; \textbf{\underline{without transient}} (matches $H_{1,1}$)'];
            figureTitle2=[figureTitle2,'; \textbf{\underline{without transient}} (matches $H_{2,1}$)'];
        end

        for ii=1:length(w_F1)
            f_rows_labels_col{ii}=['$f_{1} (t)\;,:w_{1}=',num2str(w_F1(ii)/w_d_r_vec(1)),'\ \omega_{\textrm{d},1}$'];
            w_F_col=zeros(N,1);
            w_F_col(1)=w_F1(ii);
            f_rows(ii,:)=F_0_col(1)*sin(w_F_col(1)*t_row);
            x_rows_temp=MDOF_Harmonic_Response_Visc(EigValues_mat, EigVectors_Normalized, F_0_col, w_F_col, t_row,ignoreTransient);
            x_rows1(ii,:)=x_rows_temp(1,:);
            x_rows2(ii,:)=x_rows_temp(2,:);
        end
        f_rows_labels_col{end-1}='$f_{1} (t)\;,:w_{1}=\omega_{1,1}^{\textrm{AR}}$';
        f_rows_labels_col{end}='$f_{1} (t)\;,:w_{1}=\omega_{1,2}^{\min}$';
        figure
        plotResponse_ForceVertically(t_row,x_rows1,x_rows_Y_label_col(1),f_rows,f_rows_labels_col,figureTitle1,sameScale_y1)
        figure
        plotResponse_ForceVertically(t_row,x_rows2,x_rows_Y_label_col(2),f_rows,f_rows_labels_col,figureTitle2,sameScale_y1)
    end
end

%Harmonic response 2
F_0_col=zeros(N,1);
F_0_col(2)=1;
w_F2=[w_min_12,w_AR_22,0.95*w_d_r_vec(2),w_d_r_vec(2),1.05*w_d_r_vec(2),1.5*w_d_r_vec(2)];
for ignoreTransient=ignoreTransientVector
    figureTitle1=[x_rows_Y_label_col{1},' due to $f_{2} (t)=\sin\left(w_{2}t\right)$'];
    figureTitle2=[x_rows_Y_label_col{2},' due to $f_{2} (t)=\sin\left(w_{2}t\right)$'];

    if ignoreTransient
        sameScale_y1_Vector=true;
    else
        sameScale_y1_Vector=[false,true];
    end
    for sameScale_y1=sameScale_y1_Vector
        if ignoreTransient
            figureTitle1=[figureTitle1,'; \textbf{\underline{without transient}} (matches $H_{1,2}$)'];
            figureTitle2=[figureTitle2,'; \textbf{\underline{without transient}} (matches $H_{2,2}$)'];
        end
        for ii=1:length(w_F2)
            f_rows_labels_col(ii)={['$f_{2} (t)\;,:w_{2}=',num2str(w_F2(ii)/w_d_r_vec(2)),'\ \omega_{\textrm{d},2}$']};
            w_F_col=zeros(N,1);
            w_F_col(2)=w_F1(ii);
            f_rows(ii,:)=F_0_col(2)*sin(w_F_col(2)*t_row);
            x_rows_temp=MDOF_Harmonic_Response_Visc(EigValues_mat, EigVectors_Normalized, F_0_col, w_F_col, t_row,ignoreTransient);
            x_rows1(ii,:)=x_rows_temp(1,:);
            x_rows2(ii,:)=x_rows_temp(2,:);
        end
        f_rows_labels_col(1)={'$f_{2} (t)\;,:w_{1}=\omega_{1,2}^{\min}$'};
        f_rows_labels_col(2)={'$f_{2} (t)\;,:w_{1}=\omega_{2,2}^{\textrm{AR}}$'};
        figure
        plotResponse_ForceVertically(t_row,x_rows1,x_rows_Y_label_col(1),f_rows,f_rows_labels_col,figureTitle1,sameScale_y1)
        figure
        plotResponse_ForceVertically(t_row,x_rows2,x_rows_Y_label_col(2),f_rows,f_rows_labels_col,figureTitle2,sameScale_y1)
    end
end