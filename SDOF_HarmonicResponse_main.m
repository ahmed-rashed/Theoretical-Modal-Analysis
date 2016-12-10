function SDOF_HarmonicResponse_main()
set(groot,'DefaultLineLineWidth',1);

clc
close all

t_final=60;
n_points=500;

m=1;
w_n=2;
F0=1;
t_row=linspace(0,t_final,n_points);

%% Undamped SDOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta=0;
w_0_vec=[.1,.9,1,1.1,2]*w_n;
ignoreTransient=true;
x_func=@(t_row,w_0) SDOF_Harmonic_Response_Visc_mul_m(F0, w_0, w_n, zeta, t_row,ignoreTransient)/m;
F_func=@(t_row,w_0) sin(w_0*t_row);
figure;SDOF_PlotHarmonicResponse(t_row,x_func,F_func,w_n,zeta,w_0_vec,'$f(t)=\sin(\omega_{0}t)$','$f(t)$',{'$x(t)$'},false,ignoreTransient);
figure;SDOF_PlotHarmonicResponse(t_row,x_func,F_func,w_n,zeta,w_0_vec,'$f(t)=\sin(\omega_{0}t)$','$f(t)$',{'$x(t)$'},true,ignoreTransient);

ignoreTransient=false;
x_func=@(t_row,w_0) SDOF_Harmonic_Response_Visc_mul_m(F0, w_0, w_n, zeta, t_row,ignoreTransient)/m;
figure;SDOF_PlotHarmonicResponse(t_row,x_func,F_func,w_n,zeta,w_0_vec,'$f(t)=\sin(\omega_{0}t)$','$f(t)$',{'$x(t)$'},false,ignoreTransient);
figure;SDOF_PlotHarmonicResponse(t_row,x_func,F_func,w_n,zeta,w_0_vec,'$f(t)=\sin(\omega_{0}t)$','$f(t)$',{'$x(t)$'},true,ignoreTransient);
export_figure(max(double(get(groot, 'Children')))+[-4:-1]+1,'||',{'Undamped1','Undamped2','Undamped3','Undamped4'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Damped SDOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zeta_vec=[0.1,0.3];
N_zeta=length(zeta_vec);
for n=1:N_zeta
    w_d=w_n*sqrt(1-zeta_vec(n)^2);
    w_0_vec=[.1*w_d,w_n*sqrt(1-2*zeta_vec(n)^2),w_d,1.1*w_d,2*w_d];
    ignoreTransient=false;
    x_func=@(t_row,w_0) SDOF_Harmonic_Response_Visc_mul_m(F0, w_0, w_n, zeta_vec(n), t_row,ignoreTransient)/m;
    figure;SDOF_PlotHarmonicResponse(t_row,x_func,F_func,w_n,zeta_vec(n),w_0_vec,'$f(t)=\sin(\omega_{0}t)$','$f(t)$',{'$x(t)$'},true,ignoreTransient);
    export_figure(max(double(get(groot, 'Children'))),'||',{['Damped',num2str(n)]})
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Vibration Sensor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_points=1000;
t_row=linspace(0,t_final,n_points);

Y_base=1;
sameScale_y1=false;
ignoreTransient=false;
zeta_vec=[0,0.1,1/sqrt(2)];
N_zeta=length(zeta_vec);
filenames=cell(2*N_zeta,1);
Disp_base_func=@(t_row,w_0) Y_base*sin(w_0*t_row);
Acc_base_func=@(t_row,w_0) -w_0^2*Y_base*sin(w_0*t_row);
for n=1:N_zeta
    w_d=w_n*sqrt(1-zeta_vec(n)^2);
    w_0_vec=[0.1,0.4,0.6,3,5]*w_d;

    q_func=@(t_row,w_0) w_0^2*SDOF_Harmonic_Response_Visc_mul_m(Y_base, w_0, w_n, zeta_vec(n), t_row, ignoreTransient);

    %Vibrometer response
    figure
    SDOF_PlotHarmonicResponse(t_row,@(t_row,w_0) q_func(t_row+1.5*2*pi/w_0,w_0),Disp_base_func,w_n,zeta_vec(n),w_0_vec,'$y_{\textrm{Base}}(t)=\sin(\omega_{0}t)$','$y_{\textrm{Base}}(t)$',{'$q(t+1.5T_{0})$'},sameScale_y1,ignoreTransient);
    filenames{2*n-1}=['Vibrometer',int2str(n)];

    %Accelerometer response
    figure
    SDOF_PlotHarmonicResponse(t_row,@(t_row,w_0) q_func(t_row+.5*2*pi/w_0,w_0),Acc_base_func,w_n,zeta_vec(n),w_0_vec,'$y_{\textrm{Base}}(t)=\sin(\omega_{0}t)$','$\ddot{y}_{\textrm{Base}}(t)$',{'$q(t+0.5T_{0})$'},sameScale_y1,ignoreTransient);                                 
    filenames{2*n}=['Accelerometer',int2str(n)];
end
export_figure(max(double(get(groot, 'Children')))+[-2*N_zeta:-1]+1,'||',filenames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Moving vehicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y_road=1;
sameScale_y1=true;
ignoreTransient=false;
zeta_vec=[0,0.1,1/sqrt(2)];
N_zeta=length(zeta_vec);
filenames=cell(N_zeta,1);
y_road_func=@(t_row,w_0) Y_road*sin(w_0*t_row);
for n=1:N_zeta
    y_func=@(t_row,w_0) 2*Y_road*zeta_vec(n)*w_n*SDOF_Harmonic_Response_Visc_mul_m_dot(Y_road, w_0, w_n, zeta_vec(n), t_row, ignoreTransient) ...
                          + Y_road*w_n^2*SDOF_Harmonic_Response_Visc_mul_m(Y_road, w_0, w_n, zeta_vec(n), t_row, ignoreTransient);

    w_d=w_n*sqrt(1-zeta_vec(n)^2);
    w_0_vec=[0.5,0.9,1,1.1,sqrt(2),2.5]*w_d;
    figure
    SDOF_PlotHarmonicResponse(t_row,y_func,y_road_func,w_n,zeta_vec(n),w_0_vec,'$y_{\textrm{Road}}(t)=\sin(\omega_{0}t)$','$y_{\textrm{Road}}(t)$',{'$y(t)$'},sameScale_y1,ignoreTransient);
                                 
    filenames{n}=['Vehicle',int2str(n)];
end
export_figure(max(double(get(groot, 'Children')))+[-N_zeta:-1]+1,'||',filenames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
set(groot,'DefaultLineLineWidth','remove')