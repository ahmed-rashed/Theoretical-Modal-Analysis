clearvars
clc
close all
% figure

bRaw=true;

w_n=1;
% zeta=0.4;
zeta=0.1;

f_n=w_n/2/pi;
T_n=1/f_n;

% T=2*T_n;
T=6*T_n;
K=1e2;

[D_t,f_s,D_f]=samplingParameters_T_N(T,K);

t_col=(0:K-1).'*D_t;

y_func_vec=cell(4);
y_road_cols=zeros(K,1);T_2=.5*T_n;y_road_cols(t_col<=T_2)=1;y_func_vec{1}=@(t_col,zeta) SDOF_Vehicle_Step_Response(1,w_n,zeta,t_col)-SDOF_Vehicle_Step_Response(1,w_n,zeta,t_col-T_2);
y_road_cols=[y_road_cols,zeros(K,1)];y_road_cols(1,2)=1;y_func_vec{2}=@(t_col,zeta)SDOF_Vehicle_IRF(w_n,zeta,t_col);
y_road_cols=[y_road_cols,ones(K,1)];y_func_vec{3}=@(t_col,zeta)SDOF_Vehicle_Step_Response(1,w_n,zeta,t_col);
Omega=.9*w_n;Y0=rand;y_road_cols=[y_road_cols,Y0*sin(Omega*t_col)];y_func_vec{4}=@(t_col,zeta)SDOF_Harmonic_Response_dot_Visc_mul_m(2*Y0*zeta*w_n,Omega,w_n,zeta,t_col)+SDOF_Harmonic_Response_Visc_mul_m(Y0*w_n^2,Omega,w_n,zeta,t_col);

for iii=1:size(y_road_cols,2)
    figure
    h_col=SDOF_Vehicle_IRF(w_n,zeta,t_col);
    subplot(3,1,1)
    plot(t_col/T_n,h_col,'.-');
    ylabel('$h(t)$','interpreter','latex')
    grid
    set(gca,'XTickLabel',[]);

    subplot(3,1,2)
    plot(t_col/T_n,y_road_cols(:,iii),'.-');
    ylabel('$y_{\mathrm{R}}(t)$','interpreter','latex')
    set(gca,'XTickLabel',[]);

    [y_col_approx,t_z_col]=forcedResponse(h_col,y_road_cols(:,iii),D_t,bRaw);
    y_func=y_func_vec{iii};
    y_col_exact=y_func(t_z_col,zeta);

    ax=subplot(3,1,3);
    plot_response(t_z_col,y_func_vec{iii},zeta,"$t/T_{\mathrm{n}}\qquad,:T_{\mathrm{n}}=1/f_{\mathrm{n}}=2\pi/\omega_{\mathrm{n}}="+T_n+'$','','',1/T_n,ax,'southeast');
    ylabel('$y(t)$','interpreter','latex');
    grid on

    hold on
    plot(t_z_col/T_n,y_col_approx,'.-');

    errr_col=y_col_approx-y_col_exact;
    plot(t_z_col/T_n,errr_col)

    legend(["Theortical","Numerical","error"])

    %Adjust axes limits
    if bRaw
        for ii=1:2
            pos=get(subplot(3,1,ii),'Position');
            pos(3)=pos(3)/2;
            set(subplot(3,1,ii),'Position',pos)
        end
    end
end