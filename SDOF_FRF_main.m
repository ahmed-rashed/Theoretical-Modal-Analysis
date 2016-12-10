function SDOF_FRF_main()

clc
close all

DispMagLines=0;

r=linspace(0,2.5,1e3);
zeta_vec=[0,.1,.15,.2,.3,.5,1/sqrt(2),1,2];

%Viscous SDOF
    %Single curve
    SDOF_PlotFRF(r,@SDOF_FRF_Visc_mul_k,0.15,'$r$','\zeta','H_u \ k',30,[],[],[]);
    filenames={'SDOF_1_Visc_Hu-Nyq','SDOF_1_Visc_Hu-3D','SDOF_1_Visc_Hu-Real_Imag','SDOF_1_Visc_Hu-Mag_Phase','SDOF_1_Visc_Hu-Mag_Phase_SemiLog','SDOF_1_Visc_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);


    %Multiple curves
    zzeta=linspace(0,1/sqrt(2));    %Peaks curve
    r_peaks=sqrt(1-2*zzeta.*zzeta); %Peaks curve
    H_mag_peaks=1./(1-r_peaks.*r_peaks+2*i*zzeta.*r_peaks); %Peaks curve
    SDOF_PlotFRF(r,@SDOF_FRF_Visc_mul_k,zeta_vec,'$r$','\zeta','H_u \ k',DispMagLines,r_peaks,abs(H_mag_peaks),[]);
    filenames={'SDOF_Visc_Hu-Nyq','SDOF_Visc_Hu-3D','SDOF_Visc_Hu-Real_Imag','SDOF_Visc_Hu-Mag_Phase','SDOF_Visc_Hu-Mag_Phase_SemiLog','SDOF_Visc_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    SDOF_PlotFRF(r,@(r,zeta) i*r.*SDOF_FRF_Visc_mul_k(r,zeta),zeta_vec,'$r$','\zeta','H_v \ k/\omega_{\textrm{n}}');
    filenames={'SDOF_Visc_Hv-Nyq','SDOF_Visc_Hv-3D','SDOF_Visc_Hv-Real_Imag','SDOF_Visc_Hv-Mag_Phase','SDOF_Visc_Hv-Mag_Phase_SemiLog','SDOF_Visc_Hv-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    SDOF_PlotFRF(r,@(r,zeta) -r.*r.*SDOF_FRF_Visc_mul_k(r,zeta),zeta_vec,'$r$','\zeta','H_a \ k/\omega_{\textrm{n}}^2');
    filenames={'SDOF_Visc_Ha-Nyq','SDOF_Visc_Ha-3D','SDOF_Visc_Ha-Real_Imag','SDOF_Visc_Ha-Mag_Phase','SDOF_Visc_Ha-Mag_Phase_SemiLog','SDOF_Visc_Ha-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
    close all;
    
%MovingVehicle
zzeta=linspace(0,5);    %Peaks curve
r_peaks=sqrt((sqrt(1+8*zzeta.*zzeta)-1)/4./zzeta./zzeta);   %Peaks curve
H_mag_peaks=(2*i*zzeta.*r_peaks+1)./(1-r_peaks.*r_peaks+2*i*zzeta.*r_peaks);    %Peaks curve
SDOF_PlotFRF(r,@(r,zeta) (2*i*zeta*r+1)./(1-r.*r+2*i*zeta*r),zeta_vec,'$r$','\zeta','H_u',DispMagLines,r_peaks,abs(H_mag_peaks),sqrt(2));
filenames={'SDOF_Vehcl_Hu-Nyq','SDOF_Vehcl_Hu-3D','SDOF_Vehcl_Hu-Real_Imag','SDOF_Vehcl_Hu-Mag_Phase','SDOF_Vehcl_Hu-Mag_Phase_SemiLog','SDOF_Vehcl_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
close all;
    
%Accelerometer
zzeta=linspace(0,1/sqrt(2));    %Peaks curve
r_peaks=1./sqrt(1-2*zzeta.*zzeta);  %Peaks curve
H_mag_peaks=r_peaks.*r_peaks./(1-r_peaks.*r_peaks+2*i*zzeta.*r_peaks);  %Peaks curve
SDOF_PlotFRF(r,@(r,zeta) r.*r./(1-r.*r+2*i*zeta*r),zeta_vec,'$r$','\zeta','H_u',DispMagLines,r_peaks,abs(H_mag_peaks),[]);
filenames={'SDOF_Acc_Hu-Nyq','SDOF_Acc_Hu-3D','SDOF_Acc_Hu-Real_Imag','SDOF_Acc_Hu-Mag_Phase','SDOF_Acc_Hu-Mag_Phase_SemiLog','SDOF_Acc_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
%SDOF_PlotFRF(r,@(r,zeta) r/i./(1-r.*r+2*i*zeta*r),zeta_vec,'$r$','\zeta','H_{\textrm{Vib}}',DispMagLines,[],[],[],false);
SDOF_PlotFRF(r,@(r,zeta) -1./(1-r.*r+2*i*zeta*r),zeta_vec,'$r$','\zeta','H_{\textrm{Acc}} \ \omega^2_{\textrm{n}}',DispMagLines,[],[],[]);
filenames={'SDOF_Acc_H_Acc-Nyq','SDOF_Acc_H_Acc-3D','SDOF_Acc_H_Acc-Real_Imag','SDOF_Acc_H_Acc-Mag_Phase','SDOF_Acc_H_Acc-Mag_Phase_SemiLog','SDOF_Acc_H_Acc-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
close all;
    
%Hysteretic SDOF
SDOF_PlotFRF(r,@SDOF_FRF_Hyst_mul_k,zeta_vec,'$r$','\eta','H_u \ k',DispMagLines,[],[],[]);
filenames={'SDOF_Hyst_Hu-Nyq','SDOF_Hyst_Hu-3D','SDOF_Hyst_Hu-Real_Imag','SDOF_Hyst_Hu-Mag_Phase','SDOF_Hyst_Hu-Mag_Phase_SemiLog','SDOF_Hyst_Hu-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
SDOF_PlotFRF(r,@(r,eta) i*r.*SDOF_FRF_Hyst_mul_k(r,eta),zeta_vec,'$r$','\eta','H_v \ k/\omega_{\textrm{n}}');
filenames={'SDOF_Hyst_Hv-Nyq','SDOF_Hyst_Hv-3D','SDOF_Hyst_Hv-Real_Imag','SDOF_Hyst_Hv-Mag_Phase','SDOF_Hyst_Hv-Mag_Phase_SemiLog','SDOF_Hyst_Hv-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);
SDOF_PlotFRF(r,@(r,eta) -r.*r.*SDOF_FRF_Hyst_mul_k(r,eta),zeta_vec,'$r$','\eta','H_a \ k/\omega_{\textrm{n}}^2');
filenames={'SDOF_Hyst_Ha-Nyq','SDOF_Hyst_Ha-3D','SDOF_Hyst_Ha-Real_Imag','SDOF_Hyst_Ha-Mag_Phase','SDOF_Hyst_Ha-Mag_Phase_SemiLog','SDOF_Hyst_Ha-Mag_Phase_LogLog'};export_figure(max(double(get(groot,'Children')))+[-5:0],'',filenames);