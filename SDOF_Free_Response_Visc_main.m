function SDOF_Free_Response_Visc_main()
clc
close all

set(groot,'DefaultAxesColorOrder',[0,0,1;0,0,0;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultAxesLineStyleOrder','-|--|-.')
set(groot,'DefaultLineLineWidth',1);

w_n=1;
x0=-1;
v0=0;

zeta_vec=[0,.1,.2,.4,1/sqrt(2),1,2];
legend_string={'$\zeta = 0$','$\zeta = 0.1$','$\zeta = 0.2$','$\zeta = 0.4$','$\zeta=1/\sqrt{2}$','$\zeta = 1$','$\zeta = 2$'};

t_vec=linspace(0,4*pi,500);

figure
hold on
for n=1:length(zeta_vec)
    x_vec=SDOF_Free_Response_Visc(w_n,zeta_vec(n),x0,v0,t_vec);
    plot(w_n*t_vec,x_vec)
end

title('$x(t)$ for $\omega_{n}=1$, $x_{0}=-1$ and $\dot{x}_{0}=0$','interpreter','latex');
xlabel('$\omega_{n} t$','interpreter','latex');
legend(legend_string,'interpreter','latex','Location','SouthEast');

grid on
ax=gca;
ax.XTick=0:pi:4*pi;
ax.XTickLabel={'0','\pi','2\pi','3\pi','4\pi'};
ax.XAxis.MinorTickValues=setdiff(0:pi/2:4*pi,0:pi:4*pi);
ax.XMinorGrid='on';
ax.XLim=[0,4*pi];

export_figure(gcf,'',{'SDOF_FreeResponse'})

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesLineStyleOrder','remove')
set(groot,'DefaultLineLineWidth','remove');