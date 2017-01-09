function SDOF_PlotFRF(r_vec,H_func,zeta_vec, ...
                  r_label,zeta_subtitle,H_subtitle,DispMagLines,r_peaks,H_mag_peaks,r_special,maxPhaseLag)   %Optional arguments

if nargin<4
    r_label='$r$';
else
    if isempty(r_label)
        r_label='$r$';
    end
end

if nargin<5
    zeta_subtitle='\zeta';
else
    if isempty(zeta_subtitle)
        zeta_subtitle='\zeta';
    end
end

if nargin<6
    H_subtitle='H';
else
    if isempty(H_subtitle)
        H_subtitle='H';
    end
end

if nargin<7
    DispMagLines=false;
else
    if isempty(DispMagLines)
        DispMagLines=false;
    end
end

if nargin<11
    maxPhaseLag=[];
end

set(groot,'DefaultAxesColorOrder',[0,0,1;0,0,0;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultAxesLineStyleOrder','-|--|-.')
set(groot,'DefaultLineMarkerSize',5);
set(groot,'DefaultLineLineWidth',1);

%Identify and fix the maximum axes limits
zeta_temp=min(zeta_vec(zeta_vec>0));
H_temp_vec=H_func(r_vec,zeta_temp);

Fig_Nyq=figure;
plot_FRF_Nyq(H_temp_vec,H_subtitle);
ax_Lims_Nyq=axis;clf;
axis(ax_Lims_Nyq);hold on
real_ticks=get(cla,'XTick');
imag_ticks=get(cla,'YTick');

Fig_3D=figure;
plot_FRF_3d(r_vec,H_temp_vec);
clf;
ylim(ax_Lims_Nyq(1:2));
set(cla,'YTick',real_ticks);
zlim(ax_Lims_Nyq(3:4));
set(cla,'ZTick',imag_ticks);
hold on

Fig_r_i=figure;
ax_r=subplot(2,1,1);
ylim(ax_Lims_Nyq(1:2));
set(ax_r,'YTick',real_ticks);
hold on
ax_i=subplot(2,1,2);
ylim(ax_Lims_Nyq(3:4));
set(ax_i,'YTick',imag_ticks);
hold on

Fig_Bode1=figure;

[ax_mag1]=plot_FRF_mag_phase(r_vec,H_temp_vec,[],[],[],[],[],[],maxPhaseLag);
v=axis(ax_mag1);
clf;

ax_mag1=subplot(4,1,[1:3]);hold on
ylim(ax_mag1,v(3:4));
ax_phase1=subplot(4,1,4);hold on

Fig_Bode2=figure;
ax_mag2=subplot(4,1,[1:3]);hold on
ax_phase2=subplot(4,1,4);hold on

Fig_Bode3=figure;
ax_mag3=subplot(4,1,[1:3]);hold on
ax_phase3=subplot(4,1,4);hold on

legend_str=cell(length(zeta_vec),1);
for ii=1:length(zeta_vec)
    H_vec=H_func(r_vec,zeta_vec(ii));

    H_temp_vec=H_vec;
    if zeta_vec(ii)==0
        H_temp_vec(real(H_vec)<ax_Lims_Nyq(1))=nan+1i*nan;
        H_temp_vec(real(H_vec)>ax_Lims_Nyq(2))=nan+1i*nan;
        H_temp_vec(imag(H_vec)<ax_Lims_Nyq(3))=nan+1i*nan;
        H_temp_vec(imag(H_vec)>ax_Lims_Nyq(4))=nan+1i*nan;
    end
    
    figure(Fig_3D)
    plot_FRF_3d(r_vec,H_temp_vec,r_label,H_subtitle,1,DispMagLines);
    
    plot_FRF_r_i(r_vec,H_temp_vec,ax_r,ax_i,r_label,H_subtitle);

    figure(Fig_Nyq)
    plot_FRF_Nyq(H_temp_vec,H_subtitle);

    plot_FRF_mag_phase(r_vec,H_vec,true,ax_mag1,ax_phase1,r_label,H_subtitle,DispMagLines,maxPhaseLag);
    plot_FRF_mag_phase(r_vec,H_vec,false,ax_mag2,ax_phase2,r_label,H_subtitle,DispMagLines,maxPhaseLag);
    plot_FRF_mag_phase(r_vec,H_vec,[false,false],ax_mag3,ax_phase3,r_label,H_subtitle,DispMagLines,maxPhaseLag);
    
    if zeta_vec(ii)==0
        legend_str(ii)=cellstr(['$',zeta_subtitle,'=0;\;\mathrm{misleading}$']);
    elseif zeta_vec(ii)==1/sqrt(2)
        legend_str(ii)=cellstr(['$',zeta_subtitle,'=1/\sqrt{2}$']);
    elseif zeta_vec(ii)==sqrt(2)
        legend_str(ii)=cellstr(['$',zeta_subtitle,'=\sqrt{2}$']);
    else
        legend_str(ii)=cellstr(['$',zeta_subtitle,'=',num2str(zeta_vec(ii)),'$']);
    end
end

%Peaks curve
if nargin>8 && ~isempty(r_peaks)
    H_mag_peaks(r_peaks>max(r_vec))=nan;
    %r_peaks(r_peaks>max(r_vec))=nan;
    plot(ax_mag1,r_peaks,H_mag_peaks,'-k','LineWidth',.5);
end

if nargin>9
    if ~isempty(r_special)
        xticks=get(ax_mag1,'XTick');
        set(ax_mag1,'XTick',sort(unique([xticks,r_special])));
        set(ax_phase1,'XTick',sort(unique([xticks,r_special])));
    end
end

figure(Fig_3D)
legend(legend_str,'interpreter','latex','Location','bestOutside')

figure(Fig_Nyq)
legend(legend_str,'interpreter','latex','Location','bestOutside')

legend(ax_r,legend_str,'interpreter','latex','Location','NorthEast')

legend(ax_mag1,legend_str,'interpreter','latex','Location','NorthEast')
legend(ax_mag2,legend_str,'interpreter','latex','Location','NorthEast')
legend(ax_mag3,legend_str,'interpreter','latex','Location','NorthWest')

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultAxesLineStyleOrder','remove')
set(groot,'DefaultLineMarkerSize','remove');
set(groot,'DefaultLineLineWidth','remove');