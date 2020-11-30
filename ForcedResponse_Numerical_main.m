clearvars
clc
close all
% figure

w_n=1;
% zeta=0.4;
zeta=0.1;

f_n=w_n/2/pi;
T_n=1/f_n;

% T=2*T_n;
T=6*T_n;
K=1e2;

[Delta_t,f_s,Delta_f]=samplingParameters_T_N(T,K);

t_vec=(0:K-1)*Delta_t;

h_vec=SDOF_Vehicle_IRF(w_n,zeta,t_vec);
subplot(3,1,1)
plot(t_vec/T_n,h_vec,'.-');
ylabel('$h(t)$','interpreter','latex')
grid
set(gca,'XTickLabel',[]);

y_road_vec=zeros(1,K);T_2=.5*T_n;y_road_vec(t_vec<=T_2)=1;y_func=@(t_vec,zeta) SDOF_Vehicle_Step_Response(1,w_n,zeta,t_vec)-SDOF_Vehicle_Step_Response(1,w_n,zeta,t_vec-T_2);
% y_road_vec=zeros(1,K);y_road_vec(1)=1;y_func=@(t_vec,zeta)SDOF_Vehicle_IRF(w_n,zeta,t_vec);
% y_road_vec=ones(1,K);y_func=@(t_vec,zeta)SDOF_Vehicle_Step_Response(1,w_n,zeta,t_vec);
% w_0=.9*w_n;Y0=rand;y_road_vec=Y0*sin(w_0*t_vec);y_func=@(t_vec,zeta)SDOF_Harmonic_Response_dot_Visc_mul_m(2*Y0*zeta*w_n, w_0, w_n, zeta, t_vec)+SDOF_Harmonic_Response_Visc_mul_m(Y0*w_n^2, w_0, w_n, zeta, t_vec);

subplot(3,1,2)
plot(t_vec/T_n,y_road_vec,'.-');
ylabel('$y_{\mathrm{Road}}(t)$','interpreter','latex')
set(gca,'XTickLabel',[]);

[y_vec_approx,t_z_vec]=forcedResponse(h_vec,y_road_vec,Delta_t,true);
y_vec_exact=y_func(t_z_vec,zeta);

ax=subplot(3,1,3);
plot_response(t_z_vec,y_func,zeta,'$t/T_{\mathrm{n}}\qquad,:T_{\mathrm{n}}=1/f_{\mathrm{n}}=2\pi/\omega_{\mathrm{n}}$','','',1/T_n,ax,'southeast');
ylabel('$y(t)$','interpreter','latex');
grid on

hold on
plot(t_z_vec/T_n,y_vec_approx,'.-');

errr=y_vec_approx-y_vec_exact;
plot(t_z_vec/T_n,errr)

legend({'Theortical','Numerical','error'})

%Adjust axes limits
for ii=1:2
    pos=get(subplot(3,1,ii),'Position');
    pos(3)=pos(3)/2;
    set(subplot(3,1,ii),'Position',pos)
end