clearvars
clc
close all

t_final=60;
n_points=500;

m=1;
w_n=2;
F0=1;
t_row=linspace(0,t_final,n_points);

gg=groot;

f_str="f(t)";
x_str="x(t)";
f_title_str="$"+f_str+"=F_{0}\sin(\Omega t)=F_{0}\sin(r\omega_{\mathrm{n}}t)\quad,:r\equiv\frac{\Omega}{\omega_{\mathrm{n}}}$";

%% Undamped SDOF
zeta=0;
Omega_vec=[.1,.9,1,1.1,2]*w_n;
ignoreTransient=true;
x_func=@(t_row,Omega) SDOF_Harmonic_Response_Visc_mul_m(F0,Omega,w_n,zeta,t_row,ignoreTransient)/m;
F_func=@(t_row,Omega) sin(Omega*t_row);
figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta,Omega_vec,f_title_str,f_str,x_str,false,ignoreTransient);
figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta,Omega_vec,f_title_str,f_str,x_str,true,ignoreTransient);

x_func=@(t_row,Omega) SDOF_Harmonic_Response_Visc_mul_m(F0,Omega,w_n,zeta,t_row)/m;
figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta,Omega_vec,f_title_str,f_str,x_str,false);
figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta,Omega_vec,f_title_str,f_str,x_str,true);
figs=sort([gg.Children.Number]);export_figure(figs(end+(-4:1:-1)+1),'||',"Undamped"+(1:4));

%% Damped SDOF
zeta_vec=[0.01,0.1];
N_zeta=length(zeta_vec);
for n=1:N_zeta
%     w_H_max=sqrt(1-2*zeta_vec(n)^2)*w_n;
%     Omega_vec(2)=w_H_max;
    x_func=@(t_row,Omega) SDOF_Harmonic_Response_Visc_mul_m(F0,Omega,w_n,zeta_vec(n),t_row)/m;
    figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta_vec(n),Omega_vec,f_title_str,f_str,x_str,true);
    export_figure(gcf,'||',"Damped"+n)
end

%% Vibration Sensor

n_points=1000;
t_row=linspace(0,t_final,n_points);
Omega_vec=[0.1,0.4,0.6,3,5]*w_n;

Y0=1;
sameScale_y1=false;
zeta_vec=[0.1,1/sqrt(2)];
N_zeta=length(zeta_vec);
filenames=strings(2*N_zeta,1);
Disp_base_func=@(t_row,Omega) Y0*sin(Omega*t_row);
Acc_base_func=@(t_row,Omega) -Omega^2*Y0*sin(Omega*t_row);
fig_vec=gobjects(1,2*N_zeta);
for n=1:N_zeta
    q_func=@(t_row,Omega) SDOF_Harmonic_Response_Visc_mul_m(Omega^2*Y0,Omega,w_n,zeta_vec(n),t_row);

    %Vibrometer response
    f_str="y_{\mathrm{B}}(t)";
    x_str="q(t+1.5T_{0})";
    f_title_str="$"+f_str+"=Y_{0}\sin(\Omega t)=Y_{0}\sin(r\omega_{\mathrm{n}}t)\quad,:r\equiv\frac{\Omega}{\omega_{\mathrm{n}}}$";
    fig_vec(2*n-1)=figure;
    SDOF_Plot_Harmonic_Response(t_row,@(t_row,Omega) q_func(t_row+1.5*2*pi/Omega,Omega),Disp_base_func,w_n,zeta_vec(n),Omega_vec,f_title_str,f_str,x_str,sameScale_y1);
    filenames(2*n-1)="Vibrometer"+n;

    %Accelerometer response
    f_str="\ddot{y}_{\mathrm{B}}(t)";
    x_str="q(t+0.5T_{0})";
    f_title_str="$"+f_str+"=-\Omega^{2} Y_{0}\sin(\Omega t)=-\Omega^{2} Y_{0}\sin(r\omega_{\mathrm{n}}t)\quad,:r\equiv\frac{\Omega}{\omega_{\mathrm{n}}}$";
    fig_vec(2*n)=figure;
    SDOF_Plot_Harmonic_Response(t_row,@(t_row,Omega) q_func(t_row+.5*2*pi/Omega,Omega),Acc_base_func,w_n,zeta_vec(n),Omega_vec,f_title_str,f_str,x_str,sameScale_y1,false,false);
    filenames(2*n)="Accelerometer"+n;
end
export_figure(fig_vec,'||',filenames)

%% Moving vehicle
f_str="y_{\mathrm{R}}(t)";
f_title_str="$"+f_str+"=Y_{0}\sin(\Omega t)=Y_{0}\sin(r\omega_{\mathrm{n}}t)\quad,:r\equiv\frac{\Omega}{\omega_{\mathrm{n}}}$";

Omega_vec=[0.5,0.9,1,1.1,sqrt(2),2.5]*w_n;

Y0=1;
sameScale_y1=true;
zeta_vec=[0.1,1/sqrt(2)];
N_zeta=length(zeta_vec);
filenames=strings(2*N_zeta,1);
y_road_func=@(t_row,Omega) Y0*sin(Omega*t_row);
for n=1:N_zeta
    y_Vehicle=@(t_row,Omega) SDOF_Harmonic_Response_dot_Visc_mul_m(2*Y0*zeta_vec(n)*w_n,Omega,w_n,zeta_vec(n),t_row) ...
                          + SDOF_Harmonic_Response_Visc_mul_m(Y0*w_n^2,Omega,w_n,zeta_vec(n),t_row);

    y_Acc_Vehicle=@(t_row,Omega) SDOF_Vehicle_Harmonic_Acc_Response_Visc(Y0,Omega,w_n,zeta_vec(n),t_row);

    fig_vec(2*n-1)=figure;
    x_str="y(t)";
    SDOF_Plot_Harmonic_Response(t_row,y_Vehicle,y_road_func,w_n,zeta_vec(n),Omega_vec,f_title_str,f_str,x_str,sameScale_y1);
    filenames(2*n-1)="VehicleResponse"+n;

    fig_vec(2*n)=figure;
    x_str="\ddot{y}(t)";
    SDOF_Plot_Harmonic_Response(t_row,y_Acc_Vehicle,y_road_func,w_n,zeta_vec(n),Omega_vec,f_title_str,f_str,x_str,sameScale_y1);
    filenames(2*n)="VehicleACC"+n;
end
export_figure(fig_vec,'||',filenames)