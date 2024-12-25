function SDOF_Plot_Harmonic_Response(t_row,x_func,f_func,w_n,zeta,w_0_vec,f_title_str,f_str,x_str,sameScale_y1, ...
                                     ignoreTransient,sameScale_y2)    %Optional arguments

if ~isstring(f_str),error('f_str must be string!'),end
if ~isstring(f_title_str),error('f_title_str must be string!'),end
if ~isstring(x_str),error('x_str must be string!'),end

set(groot,'DefaultLineLineWidth',1);

if nargin<11
    ignoreTransient=false;
end

if nargin<12
    sameScale_y2=true;
end

n_points=length(t_row);
ii_row=length(w_0_vec);

if zeta==1/sqrt(2)
    zeta_expr="$\zeta=1/\sqrt{2}$";
elseif zeta==sqrt(2)
    zeta_expr="$\zeta=\sqrt{2}$";
else
    zeta_expr="$\zeta="+zeta+'$';
end

if ignoreTransient
    x_modified_str=strrep(x_str,'(','_{\mathrm{ss}}(');
else
    x_modified_str=x_str;
end

x_title_str="$"+x_modified_str+'$';
figureTitle=zeta_expr;
if ignoreTransient && zeta==0
    figureTitle=zeta_expr+". "+x_title_str+' never coincides with $'+x_str+'$';
end

r_str_col=strings(length(ii_row),1);
x_rows=zeros(ii_row,n_points);
f_rows=zeros(ii_row,n_points);
for ii=1:ii_row
    x_rows(ii,:)=x_func(t_row,w_0_vec(ii));
    f_rows(ii,:)=f_func(t_row,w_0_vec(ii));
    r_str_col(ii)="$r="+(w_0_vec(ii)/w_n)+'$';
end

plot_Forced_Response_Vertically(t_row,x_rows,f_rows,figureTitle,f_title_str,r_str_col,x_title_str,sameScale_y1,sameScale_y2);

set(groot,'DefaultLineLineWidth','remove')