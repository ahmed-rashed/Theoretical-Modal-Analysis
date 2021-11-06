clearvars
clc
close all

DispMagLines=0;

r_vec=linspace(0,2.5,1e3);
zeta_vec=[0,.1,.15,.2,.3,.5,1/sqrt(2),1,2];

%% Viscous SDOF
%Single curve
SDOF_PlotFRF(r_vec,@SDOF_FRF_Visc_mul_k,0.15,'','',["H_{u}","k"],30);
filenames=["SDOF_1_Visc_Hu-Nyq","SDOF_1_Visc_Hu-3D","SDOF_1_Visc_Hu-Real_Imag","SDOF_1_Visc_Hu-Mag_Phase","SDOF_1_Visc_Hu-Mag_Phase_SemiLog","SDOF_1_Visc_Hu-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);


%Multiple curves
zzeta=linspace(0,1/sqrt(2));    %Peaks curve
r_peaks=sqrt(1-2*zzeta.*zzeta); %Peaks curve
H_mag_peaks=1./(1-r_peaks.*r_peaks+2*1i*zzeta.*r_peaks); %Peaks curve
SDOF_PlotFRF(r_vec,@SDOF_FRF_Visc_mul_k,zeta_vec,'','',["H_{u}","k"],DispMagLines,r_peaks,abs(H_mag_peaks));
filenames=["SDOF_Visc_Hu-Nyq","SDOF_Visc_Hu-3D","SDOF_Visc_Hu-Real_Imag","SDOF_Visc_Hu-Mag_Phase","SDOF_Visc_Hu-Mag_Phase_SemiLog","SDOF_Visc_Hu-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
SDOF_PlotFRF(r_vec,@(r_vec,zeta) 1i*r_vec.*SDOF_FRF_Visc_mul_k(r_vec,zeta),zeta_vec,'','',["H_{v}","k/\omega_{\mathrm{n}}"]);
filenames=["SDOF_Visc_Hv-Nyq","SDOF_Visc_Hv-3D","SDOF_Visc_Hv-Real_Imag","SDOF_Visc_Hv-Mag_Phase","SDOF_Visc_Hv-Mag_Phase_SemiLog","SDOF_Visc_Hv-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
SDOF_PlotFRF(r_vec,@(r_vec,zeta) -r_vec.^2.*SDOF_FRF_Visc_mul_k(r_vec,zeta),zeta_vec,'','',["H_{a}","k/\omega_{\mathrm{n}}^2"]);
filenames=["SDOF_Visc_Ha-Nyq","SDOF_Visc_Ha-3D","SDOF_Visc_Ha-Real_Imag","SDOF_Visc_Ha-Mag_Phase","SDOF_Visc_Ha-Mag_Phase_SemiLog","SDOF_Visc_Ha-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    
%% MovingVehicle
r_vec=linspace(0,4,1e3);

close all;
zzeta=linspace(0,5);    %Peaks curve
r_peaks=sqrt((sqrt(1+8*zzeta.*zzeta)-1)/4./zzeta./zzeta);   %Peaks curve
H_mag_peaks=(2*1i*zzeta.*r_peaks+1)./(1-r_peaks.*r_peaks+2*1i*zzeta.*r_peaks);    %Peaks curve
fn_SDOF_vehicle_FRF=@(r_vec,zeta) (2*1i*zeta*r_vec+1)./(1-r_vec.*r_vec+2*1i*zeta*r_vec);

SDOF_PlotFRF(r_vec,fn_SDOF_vehicle_FRF,zeta_vec,'','','H_{u}',DispMagLines,r_peaks,abs(H_mag_peaks));
filenames=["SDOF_Vehcl_Hu-Nyq","SDOF_Vehcl_Hu-3D","SDOF_Vehcl_Hu-Real_Imag","SDOF_Vehcl_Hu-Mag_Phase","SDOF_Vehcl_Hu-Mag_Phase_SemiLog","SDOF_Vehcl_Hu-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);

SDOF_PlotFRF(r_vec,@(r_vec,zeta) 1i*r_vec.*fn_SDOF_vehicle_FRF(r_vec,zeta),zeta_vec,'','','H_{v}',[],[]);
filenames=["SDOF_Vehcl_Hv-Nyq","SDOF_Vehcl_Hv-3D","SDOF_Vehcl_Hv-Real_Imag","SDOF_Vehcl_Hv-Mag_Phase","SDOF_Vehcl_Hv-Mag_Phase_SemiLog","SDOF_Vehcl_Hv-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);

SDOF_PlotFRF(r_vec,@(r_vec,zeta) -r_vec.^2.*fn_SDOF_vehicle_FRF(r_vec,zeta),zeta_vec,'','','H_{a}',[],[]);
filenames=["SDOF_Vehcl_Ha-Nyq","SDOF_Vehcl_Ha-3D","SDOF_Vehcl_Ha-Real_Imag","SDOF_Vehcl_Ha-Mag_Phase","SDOF_Vehcl_Ha-Mag_Phase_SemiLog","SDOF_Vehcl_Ha-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);

%% Accelerometer
r_vec=linspace(0,2.5,1e3);

close all;
zzeta=linspace(0,1/sqrt(2));    %Peaks curve
r_peaks=1./sqrt(1-2*zzeta.*zzeta);  %Peaks curve
H_mag_peaks=r_peaks.^2./(1-r_peaks.*r_peaks+2*1i*zzeta.*r_peaks);  %Peaks curve
fn_SDOF_sensor_FRF=@(r_vec,zeta) r_vec.^2./(1-r_vec.*r_vec+2*1i*zeta*r_vec);

SDOF_PlotFRF(r_vec,fn_SDOF_sensor_FRF,zeta_vec                                      ,'','','H_{u}',DispMagLines,r_peaks,abs(H_mag_peaks));
filenames=["SDOF_Sensor_Hu-Nyq","SDOF_Sensor_Hu-3D","SDOF_Sensor_Hu-Real_Imag","SDOF_Sensor_Hu-Mag_Phase","SDOF_Sensor_Hu-Mag_Phase_SemiLog","SDOF_Sensor_Hu-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);

SDOF_PlotFRF(r_vec,@(r_vec,zeta) fn_SDOF_sensor_FRF(r_vec,zeta)/1i./r_vec,zeta_vec  ,'','',["H_{\mathrm{Vel}}","\omega_{\mathrm{n}}"]       ,DispMagLines);
filenames=["SDOF_Sensor_H_Vel-Nyq","SDOF_Sensor_H_Vel-3D","SDOF_Sensor_H_Vel-Real_Imag","SDOF_Sensor_H_Vel-Mag_Phase","SDOF_Sensor_H_Vel-Mag_Phase_SemiLog","SDOF_Sensor_H_Vel-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);

SDOF_PlotFRF(r_vec,@(r_vec,zeta) fn_SDOF_sensor_FRF(r_vec,zeta)./-r_vec.^2,zeta_vec               ,'','',["H_{\mathrm{Acc}}","\omega^{2}_{\mathrm{n}}"]   ,DispMagLines);
filenames=["SDOF_Sensor_H_Acc-Nyq","SDOF_Sensor_H_Acc-3D","SDOF_Sensor_H_Acc-Real_Imag","SDOF_Sensor_H_Acc-Mag_Phase","SDOF_Sensor_H_Acc-Mag_Phase_SemiLog","SDOF_Sensor_H_Acc-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    
%% Structural SDOF
close all;
SDOF_PlotFRF(r_vec,@SDOF_FRF_Struc_mul_k,zeta_vec,'','\eta',["H_{u}","k"],DispMagLines);
filenames=["SDOF_Struc_Hu-Nyq","SDOF_Struc_Hu-3D","SDOF_Struc_Hu-Real_Imag","SDOF_Struc_Hu-Mag_Phase","SDOF_Struc_Hu-Mag_Phase_SemiLog","SDOF_Struc_Hu-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
SDOF_PlotFRF(r_vec,@(r_vec,eta) 1i*r_vec.*SDOF_FRF_Struc_mul_k(r_vec,eta),zeta_vec,'','\eta',["H_v","k/\omega_{\mathrm{n}}"]);
filenames=["SDOF_Struc_Hv-Nyq","SDOF_Struc_Hv-3D","SDOF_Struc_Hv-Real_Imag","SDOF_Struc_Hv-Mag_Phase","SDOF_Struc_Hv-Mag_Phase_SemiLog","SDOF_Struc_Hv-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
SDOF_PlotFRF(r_vec,@(r_vec,eta) -r_vec.*r_vec.*SDOF_FRF_Struc_mul_k(r_vec,eta),zeta_vec,'','\eta',["H_a","k/\omega_{\mathrm{n}}^2"]);
filenames=["SDOF_Struc_Ha-Nyq","SDOF_Struc_Ha-3D","SDOF_Struc_Ha-Real_Imag","SDOF_Struc_Ha-Mag_Phase","SDOF_Struc_Ha-Mag_Phase_SemiLog","SDOF_Struc_Ha-Mag_Phase_LogLog"];export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
