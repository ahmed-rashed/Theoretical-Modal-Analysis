function SDOF_Plot_Harmonic_Response(t_row,x_func,f_func,w_n,zeta,w_0_vec,HarmonicExcitation_title,f_label,x_label_rows,sameScale_y1,ignoreTransient)
set(groot,'DefaultLineLineWidth',1);

n_points=length(t_row);
ii_row=length(w_0_vec);

if zeta==1/sqrt(2)
    zeta_expr='$\zeta=1/\sqrt{2}$';
elseif zeta==sqrt(2)
    zeta_expr='$\zeta=\sqrt{2}$';
else
    zeta_expr=['$\zeta=',num2str(zeta),'$'];
end

if ignoreTransient
    x_new_label_rows=strrep(x_label_rows,'(','_{\mathrm{ss}}(');
else
    x_new_label_rows=x_label_rows;
end

figureTitle=['$',x_new_label_rows,'$ due to $',HarmonicExcitation_title,'$ for ',zeta_expr];
if ignoreTransient && zeta==0
    figureTitle=[figureTitle,' \underline{(never coincides with $',x_label_rows,'$, but matches $H(\omega)$)}'];
end

f_rows_labels_col=cell(length(ii_row),1);
x_rows=zeros(ii_row,n_points);
f_rows=zeros(ii_row,n_points);
for ii=1:ii_row
    x_rows(ii,:)=x_func(t_row,w_0_vec(ii));
    f_rows(ii,:)=f_func(t_row,w_0_vec(ii));
    f_rows_labels_col(ii)={['$',f_label,',:r_{0}=',num2str(w_0_vec(ii)/w_n),'$']};
end

plot_Forced_Response_Vertically(t_row,x_rows,x_new_label_rows,f_rows,f_rows_labels_col,figureTitle,sameScale_y1);

set(groot,'DefaultLineLineWidth','remove')