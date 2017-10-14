clearvars
clc
close all

DispMagLines=0;

r_vec=linspace(0,2.5,1e3);
zeta_vec=[0,.1,.15,.2,.3,.5,1/sqrt(2),1,2];

%% Viscous SDOF
    %Single curve
    SDOF_PlotFRF(r_vec,SDOF_FRF_Visc_mul_k,0.15,'','',{'H_{u}','k'},30);
    filenames={'SDOF_1_Visc_Hu-Nyq','SDOF_1_Visc_Hu-3D','SDOF_1_Visc_Hu-Real_Imag','SDOF_1_Visc_Hu-Mag_Phase','SDOF_1_Visc_Hu-Mag_Phase_SemiLog','SDOF_1_Visc_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);


    %Multiple curves
    zzeta=linspace(0,1/sqrt(2));    %Peaks curve
    r_peaks=sqrt(1-2*zzeta.*zzeta); %Peaks curve
    H_mag_peaks=1./(1-r_peaks.*r_peaks+2*1i*zzeta.*r_peaks); %Peaks curve
    SDOF_PlotFRF(r_vec,@SDOF_FRF_Visc_mul_k,zeta_vec,'','',{'H_{u}','k'},DispMagLines,r_peaks,abs(H_mag_peaks));
    filenames={'SDOF_Visc_Hu-Nyq','SDOF_Visc_Hu-3D','SDOF_Visc_Hu-Real_Imag','SDOF_Visc_Hu-Mag_Phase','SDOF_Visc_Hu-Mag_Phase_SemiLog','SDOF_Visc_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    SDOF_PlotFRF(r_vec,@(r_vec,zeta) 1i*r_vec.*SDOF_FRF_Visc_mul_k(r_vec,zeta),zeta_vec,'','',{'H_v','k/\omega_{\mathrm{n}}'});
    filenames={'SDOF_Visc_Hv-Nyq','SDOF_Visc_Hv-3D','SDOF_Visc_Hv-Real_Imag','SDOF_Visc_Hv-Mag_Phase','SDOF_Visc_Hv-Mag_Phase_SemiLog','SDOF_Visc_Hv-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    SDOF_PlotFRF(r_vec,@(r_vec,zeta) -r_vec.*r_vec.*SDOF_FRF_Visc_mul_k(r_vec,zeta),zeta_vec,'','',{'H_a','k/\omega_{\mathrm{n}}^2'});
    filenames={'SDOF_Visc_Ha-Nyq','SDOF_Visc_Ha-3D','SDOF_Visc_Ha-Real_Imag','SDOF_Visc_Ha-Mag_Phase','SDOF_Visc_Ha-Mag_Phase_SemiLog','SDOF_Visc_Ha-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    
%% MovingVehicle
close all;
zzeta=linspace(0,5);    %Peaks curve
r_peaks=sqrt((sqrt(1+8*zzeta.*zzeta)-1)/4./zzeta./zzeta);   %Peaks curve
H_mag_peaks=(2*1i*zzeta.*r_peaks+1)./(1-r_peaks.*r_peaks+2*1i*zzeta.*r_peaks);    %Peaks curve
SDOF_PlotFRF(r_vec,@(r_vec,zeta) (2*1i*zeta*r_vec+1)./(1-r_vec.*r_vec+2*1i*zeta*r_vec),zeta_vec,'','','',DispMagLines,r_peaks,abs(H_mag_peaks));
filenames={'SDOF_Vehcl_Hu-Nyq','SDOF_Vehcl_Hu-3D','SDOF_Vehcl_Hu-Real_Imag','SDOF_Vehcl_Hu-Mag_Phase','SDOF_Vehcl_Hu-Mag_Phase_SemiLog','SDOF_Vehcl_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    
%% Accelerometer
close all;
zzeta=linspace(0,1/sqrt(2));    %Peaks curve
r_peaks=1./sqrt(1-2*zzeta.*zzeta);  %Peaks curve
H_mag_peaks=r_peaks.*r_peaks./(1-r_peaks.*r_peaks+2*1i*zzeta.*r_peaks);  %Peaks curve
SDOF_PlotFRF(r_vec,@(r_vec,zeta) r_vec.*r_vec./(1-r_vec.*r_vec+2*1i*zeta*r_vec),zeta_vec,'','','H_{u}',DispMagLines,r_peaks,abs(H_mag_peaks));
filenames={'SDOF_Acc_Hu-Nyq','SDOF_Acc_Hu-3D','SDOF_Acc_Hu-Real_Imag','SDOF_Acc_Hu-Mag_Phase','SDOF_Acc_Hu-Mag_Phase_SemiLog','SDOF_Acc_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
%SDOF_PlotFRF(r_vec,@(r_vec,zeta) r_vec/1i./(1-r_vec.*r_vec+2*1i*zeta*r_vec),zeta_vec,'','','H_{\mathrm{Vib}}',DispMagLines,[],[],[],false);
SDOF_PlotFRF(r_vec,@(r_vec,zeta) -1./(1-r_vec.*r_vec+2*1i*zeta*r_vec),zeta_vec,'','',{'H_{\mathrm{Acc}}','\omega^2_{\mathrm{n}}'},DispMagLines);
filenames={'SDOF_Acc_H_Acc-Nyq','SDOF_Acc_H_Acc-3D','SDOF_Acc_H_Acc-Real_Imag','SDOF_Acc_H_Acc-Mag_Phase','SDOF_Acc_H_Acc-Mag_Phase_SemiLog','SDOF_Acc_H_Acc-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    
%% Hysteretic SDOF
close all;
SDOF_PlotFRF(r_vec,@SDOF_FRF_Hyst_mul_k,zeta_vec,'','\eta',{'H_{u}','k'},DispMagLines);
filenames={'SDOF_Hyst_Hu-Nyq','SDOF_Hyst_Hu-3D','SDOF_Hyst_Hu-Real_Imag','SDOF_Hyst_Hu-Mag_Phase','SDOF_Hyst_Hu-Mag_Phase_SemiLog','SDOF_Hyst_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
SDOF_PlotFRF(r_vec,@(r_vec,eta) 1i*r_vec.*SDOF_FRF_Hyst_mul_k(r_vec,eta),zeta_vec,'','\eta',{'H_v','k/\omega_{\mathrm{n}}'});
filenames={'SDOF_Hyst_Hv-Nyq','SDOF_Hyst_Hv-3D','SDOF_Hyst_Hv-Real_Imag','SDOF_Hyst_Hv-Mag_Phase','SDOF_Hyst_Hv-Mag_Phase_SemiLog','SDOF_Hyst_Hv-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
SDOF_PlotFRF(r_vec,@(r_vec,eta) -r_vec.*r_vec.*SDOF_FRF_Hyst_mul_k(r_vec,eta),zeta_vec,'','\eta',{'H_a','k/\omega_{\mathrm{n}}^2'});
filenames={'SDOF_Hyst_Ha-Nyq','SDOF_Hyst_Ha-3D','SDOF_Hyst_Ha-Real_Imag','SDOF_Hyst_Ha-Mag_Phase','SDOF_Hyst_Ha-Mag_Phase_SemiLog','SDOF_Hyst_Ha-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);