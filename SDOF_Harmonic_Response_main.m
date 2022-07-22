clearvars
clc
close all

t_final=60;
n_points=500;

m=1;
w_n=2;
F0=1;
t_row=linspace(0,t_final,n_points);

%% Undamped SDOF
zeta=0;
w_0_vec=[.1,.9,1,1.1,2]*w_n;
ignoreTransient=true;
x_func=@(t_row,w_0) SDOF_Harmonic_Response_Visc_mul_m(F0,w_0,w_n,zeta,t_row,ignoreTransient)/m;
F_func=@(t_row,w_0) sin(w_0*t_row);
figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta,w_0_vec,'f(t)=\sin(\omega_{0}t)','f(t)','x(t)',false,ignoreTransient);
figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta,w_0_vec,'f(t)=\sin(\omega_{0}t)','f(t)','x(t)',true,ignoreTransient);

x_func=@(t_row,w_0) SDOF_Harmonic_Response_Visc_mul_m(F0,w_0,w_n,zeta,t_row)/m;
figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta,w_0_vec,'f(t)=\sin(\omega_{0}t)','f(t)','x(t)',false);
figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta,w_0_vec,'f(t)=\sin(\omega_{0}t)','f(t)','x(t)',true);
export_figure(max(double(get(groot,'Children')))+(-4:-1)+1,'||',"Undamped"+(1:4))

%% Damped SDOF
zeta_vec=[0.01,0.1];
N_zeta=length(zeta_vec);
for n=1:N_zeta
%     w_H_max=sqrt(1-2*zeta_vec(n)^2)*w_n;
%     w_0_vec(2)=w_H_max;
    x_func=@(t_row,w_0) SDOF_Harmonic_Response_Visc_mul_m(F0,w_0,w_n,zeta_vec(n),t_row)/m;
    figure;SDOF_Plot_Harmonic_Response(t_row,x_func,F_func,w_n,zeta_vec(n),w_0_vec,'f(t)=\sin(\omega_{0}t)','f(t)','x(t)',true);
    export_figure(max(double(get(groot,'Children'))),'||',"Damped"+n)
end

%% Vibration Sensor
n_points=1000;
t_row=linspace(0,t_final,n_points);
w_0_vec=[0.1,0.4,0.6,3,5]*w_n;

Y0=1;
sameScale_y1=false;
zeta_vec=[0.1,1/sqrt(2)];
N_zeta=length(zeta_vec);
filenames=strings(2*N_zeta,1);
Disp_base_func=@(t_row,w_0) Y0*sin(w_0*t_row);
Acc_base_func=@(t_row,w_0) -w_0^2*Y0*sin(w_0*t_row);
for n=1:N_zeta
    q_func=@(t_row,w_0) SDOF_Harmonic_Response_Visc_mul_m(w_0^2*Y0,w_0,w_n,zeta_vec(n),t_row);

    %Vibrometer response
    figure
    SDOF_Plot_Harmonic_Response(t_row,@(t_row,w_0) q_func(t_row+1.5*2*pi/w_0,w_0),Disp_base_func,w_n,zeta_vec(n),w_0_vec,'y_{\mathrm{B}}(t)=\sin(\omega_{0}t)','y_{\mathrm{B}}(t)','q(t+1.5T_{0})',sameScale_y1);
    filenames(2*n-1)="Vibrometer"+n;

    %Accelerometer response
    figure
    SDOF_Plot_Harmonic_Response(t_row,@(t_row,w_0) q_func(t_row+.5*2*pi/w_0,w_0),Acc_base_func,w_n,zeta_vec(n),w_0_vec,'y_{\mathrm{B}}(t)=\sin(\omega_{0}t)','\ddot{y}_{\mathrm{B}}(t)','q(t+0.5T_{0})',sameScale_y1);                                 
    filenames(2*n)="Accelerometer"+n;
end
export_figure(max(double(get(groot,'Children')))+(-2*N_zeta:-1)+1,'||',filenames)

%% Moving vehicle
w_0_vec=[0.5,0.9,1,1.1,sqrt(2),2.5]*w_n;

Y0=1;
sameScale_y1=true;
zeta_vec=[0.1,1/sqrt(2)];
N_zeta=length(zeta_vec);
filenames=strings(2*N_zeta,1);
y_road_func=@(t_row,w_0) Y0*sin(w_0*t_row);
for n=1:N_zeta
    y_Vehicle=@(t_row,w_0) SDOF_Harmonic_Response_dot_Visc_mul_m(2*Y0*zeta_vec(n)*w_n,w_0,w_n,zeta_vec(n),t_row) ...
                          + SDOF_Harmonic_Response_Visc_mul_m(Y0*w_n^2,w_0,w_n,zeta_vec(n),t_row);

    y_Acc_Vehicle=@(t_row,w_0) SDOF_Vehicle_Harmonic_Acc_Response_Visc(Y0,w_0,w_n,zeta_vec(n),t_row);

    figure
    SDOF_Plot_Harmonic_Response(t_row,y_Vehicle,y_road_func,w_n,zeta_vec(n),w_0_vec,'y_{\mathrm{R}}(t)=\sin(\omega_{0}t)','y_{\mathrm{R}}(t)','y(t)',sameScale_y1);
    filenames(2*n-1)="VehicleResponse"+n;

    figure
    SDOF_Plot_Harmonic_Response(t_row,y_Acc_Vehicle,y_road_func,w_n,zeta_vec(n),w_0_vec,'y_{\mathrm{R}}(t)=\sin(\omega_{0}t)','y_{\mathrm{R}}(t)','\ddot{y}(t)',sameScale_y1);
    filenames(2*n)="VehicleACC"+n;
end
export_figure(max(double(get(groot,'Children')))+(-2*N_zeta:-1)+1,'||',filenames)