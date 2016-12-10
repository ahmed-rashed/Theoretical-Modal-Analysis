function SDOF_PlotHarmonicResponse(t_row,x_func,f_func,w_n,zeta,w_0_vec,HarmonicExcitation_label,f_label,x_rows_label,sameScale_y1,ignoreTransient)

n_points=length(t_row);
ii_row=length(w_0_vec);
w_d=w_n*sqrt(1-zeta^2);

figureTitle=[x_rows_label{1},' due to $',HarmonicExcitation_label,'$ for $\zeta=',num2str(zeta),'$'];
if ignoreTransient
    figureTitle=[figureTitle,'; \textbf{\underline{without transient}} (matches $H(\omega)$)'];
end

f_rows_labels_col=cell(length(ii_row),1);
x_rows=zeros(ii_row,n_points);
f_rows=zeros(ii_row,n_points);
for ii=1:ii_row
    x_rows(ii,:)=x_func(t_row,w_0_vec(ii));
    f_rows(ii,:)=f_func(t_row,w_0_vec(ii));
    if zeta == 0
        f_rows_labels_col(ii)={[f_label,' $,:\frac{\omega_{0}}{\omega_{\textrm{n}}}=',num2str(w_0_vec(ii)/w_d),'$']};
    else
        f_rows_labels_col(ii)={[f_label,' $,:\frac{\omega_{0}}{\omega_{\textrm{d}}}=',num2str(w_0_vec(ii)/w_d),'$']};
    end
end

plotResponse_ForceVertically(t_row,x_rows,x_rows_label,f_rows,f_rows_labels_col,figureTitle,sameScale_y1);