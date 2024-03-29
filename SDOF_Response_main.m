clearvars
clc
close all

x0=-1;
x_dot_0=0;
w_n=1;
f_n=w_n/2/pi;
T_n=1/f_n;
t_vec=linspace(0,2*T_n,500);
zeta_vec=[0,.1,.2,.4,1/sqrt(2),1,2];

%% Free response of viscous SDOF
plot_response(t_vec,@(t_vec,zeta)SDOF_Free_Response_Visc(w_n,zeta,x0,x_dot_0,t_vec),zeta_vec,'$t/T_{\mathrm{n}}\qquad,:T_{\mathrm{n}}=1/f_{\mathrm{n}}=2\pi/\omega_{\mathrm{n}}$','','',1/T_n,[],'southeast');
title("$x(t)$ for $\omega_{n}="+w_n+'$, $x_{0}='+x0+'$ and $\dot{x}_{0}='+x_dot_0+'$','interpreter','latex');
grid on

export_figure(gcf,'',"SDOF_FreeResponse")

%% Viscous SDOF IRF
plot_response(t_vec,@(t_vec,zeta)SDOF_IRF_Visc_mul_m(w_n,zeta,t_vec),zeta_vec,'$t/T_{\mathrm{n}}\qquad,:T_{\mathrm{n}}=1/f_{\mathrm{n}}=2\pi/\omega_{\mathrm{n}}$','','',1/T_n);
title("$h(t) \ m$ for $\omega_{n}="+w_n+'$','interpreter','latex');
grid on
export_figure(gcf,'',"SDOF_IRF")

%% Moving Vehicle IRF
plot_response(t_vec,@(t_vec,zeta)SDOF_Vehicle_IRF(w_n,zeta,t_vec),zeta_vec,'$t/T_{\mathrm{n}}\qquad,:T_{\mathrm{n}}=1/f_{\mathrm{n}}=2\pi/\omega_{\mathrm{n}}$','','',1/T_n);
title("Vehicle''s $h(t)$ for $\omega_{n}="+w_n+'$','interpreter','latex');
grid on
export_figure(gcf,'',"Vehicle_IRF")

%% Moving Vehicle Step Response
plot_response(t_vec,@(t_vec,zeta)SDOF_Vehicle_Step_Response(1,w_n,zeta,t_vec),zeta_vec,'$t/T_{\mathrm{n}}\qquad,:T_{\mathrm{n}}=1/f_{\mathrm{n}}=2\pi/\omega_{\mathrm{n}}$','','',1/T_n,[],'southeast');
title("$y_{\mathrm{step}}(t)$ for $\omega_{n}="+w_n+'$','interpreter','latex');
grid on
export_figure(gcf,'',"Vehicle_StepResponse")

%% Moving Vehicle Obstacle Response
T_2=0.05*T_n;
plot_response(t_vec,@(t_vec,zeta) SDOF_Vehicle_Step_Response(1,w_n,zeta,t_vec)-SDOF_Vehicle_Step_Response(1,w_n,zeta,t_vec-T_2),zeta_vec,'$t/T_{\mathrm{n}}\qquad,:T_{\mathrm{n}}=1/f_{\mathrm{n}}=2\pi/\omega_{\mathrm{n}}$','','',1/T_n,[],'northeast');
title("$y_{\mathrm{step}}\left(t\right)-y_{\mathrm{step}}\left(t-"+(T_2/T_n)+'T_{\mathrm{n}}\right)$ for $\omega_{n}='+w_n+'$','interpreter','latex');
grid on
% set(gca,'XAxisLocation','origin')
ylim(2*[-1,1])
export_figure(gcf,'',"Vehicle_StepResponse_1")

T_2=.5*T_n;
plot_response(t_vec,@(t_vec,zeta) SDOF_Vehicle_Step_Response(1,w_n,zeta,t_vec)-SDOF_Vehicle_Step_Response(1,w_n,zeta,t_vec-T_2),zeta_vec,'$t/T_{\mathrm{n}}\qquad,:T_{\mathrm{n}}=1/f_{\mathrm{n}}=2\pi/\omega_{\mathrm{n}}$','','',1/T_n,[],'southeast');
title("$y_{\mathrm{step}}\left(t\right)-y_{\mathrm{step}}\left(t-"+(T_2/T_n)+'T_{\mathrm{n}}\right)$ for $\omega_{n}='+w_n+'$','interpreter','latex');
grid on
% set(gca,'XAxisLocation','origin')
export_figure(gcf,'',"Vehicle_StepResponse_2")