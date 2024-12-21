function SDOF_PlotFRF(r_vec,H_func,zeta_vec, ...
                  r_label,zeta_subtitle,H_subtitle,DispMagLines,r_peaks,H_mag_peaks,r_special,maxPhaseLag)   %Optional arguments

if nargin<4
    r_label='$r\equiv\frac{\omega}{\omega_{\mathrm{n}}}$';
else
    if isempty(r_label)
        r_label='$r\equiv\frac{\omega}{\omega_{\mathrm{n}}}$';
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

set(groot,'DefaultAxesLineStyleOrder','-|--|-.')
set(groot,'DefaultLineMarkerSize',5);

N_zeta=length(zeta_vec);
Fig_Nyq=figure;
Fig_3D=figure;
Fig_r_i=figure;
Fig_mag1=figure;
if N_zeta>1
    %Identify and fix the maximum axes limits
    zeta_temp=min(zeta_vec(zeta_vec>0));
    H_temp_vec=H_func(r_vec,zeta_temp);
    
    figure(Fig_Nyq);
    plot_FRF_Nyq(H_temp_vec);
    axis("padded")
    ax_Lims_Nyq=axis;
    real_ticks=get(gca,'XTick');
    imag_ticks=get(gca,'YTick');
    
    figure(Fig_mag1); % bode plot
    ax_mag1=plot_FRF_mag_phase(r_vec,H_temp_vec,[],[],[],[],[],[],maxPhaseLag);
    v=axis(ax_mag1);
    r_ticks=get(ax_mag1,'XTick');
    clf

    figure(Fig_3D)
    plot_FRF_3d(r_vec,H_temp_vec);
    % Fig_3D_AR=get(gca,'DataAspectRatio');
    xlim(v(1:2));
    ylim(ax_Lims_Nyq(1:2));
    zlim(ax_Lims_Nyq(3:4));
    set(gca,'XTick',r_ticks,'YTick',real_ticks,'ZTick',imag_ticks);
    cla
    hold on
end

ax_r=[];
ax_i=[];
ax_mag1=[];
ax_phase1=[];
ax_mag2=[];
ax_phase2=[];
ax_mag3=[];
ax_phase3=[];
legend_str_vec=strings(N_zeta,1);
for ii=1:N_zeta
    H_vec=H_func(r_vec,zeta_vec(ii));

    H_temp_vec=H_vec;
    if N_zeta>1 && zeta_vec(ii)==0
        [~,ind]=max(abs(H_vec));
        H_temp_vec(ind)=nan+1i*nan;
    end

    figure(Fig_3D)
    plot_FRF_3d(r_vec,H_temp_vec,r_label,H_subtitle,1,DispMagLines);

    if ii==1,figure(Fig_r_i);end
    [ax_r,ax_i]=plot_FRF_r_i(r_vec,H_temp_vec,ax_r,ax_i,r_label,H_subtitle);
    
    figure(Fig_Nyq)
    plot_FRF_Nyq(H_temp_vec,r_label,H_subtitle);

    if ii==1,figure(Fig_mag1),end
    [ax_mag1,ax_phase1]=plot_FRF_mag_phase(r_vec,H_vec,true,ax_mag1,ax_phase1,r_label,H_subtitle,DispMagLines,maxPhaseLag);
    if ii==1,figure;end
    [ax_mag2,ax_phase2]=plot_FRF_mag_phase(r_vec,H_vec,false,ax_mag2,ax_phase2,r_label,H_subtitle,DispMagLines,maxPhaseLag);
    if ii==1,figure;end
    [ax_mag3,ax_phase3]=plot_FRF_mag_phase(r_vec,H_vec,[false,false],ax_mag3,ax_phase3,r_label,H_subtitle,DispMagLines,maxPhaseLag);

    if ii==1
        figure(Fig_3D)
        hold on
        figure(Fig_Nyq);
        hold on
        hold([ax_r,ax_i],"on");
        hold([ax_mag1,ax_phase1],"on");
        hold([ax_mag2,ax_phase2],"on");
        hold([ax_mag3,ax_phase3],"on");
    end

    if zeta_vec(ii)==0
        legend_str_vec(ii)="$"+zeta_subtitle+'=0;\;\mathrm{misleading}$';
    elseif zeta_vec(ii)==1/sqrt(2)
        legend_str_vec(ii)="$"+zeta_subtitle+'=1/\sqrt{2}$';
    elseif zeta_vec(ii)==sqrt(2)
        legend_str_vec(ii)="$"+zeta_subtitle+'=\sqrt{2}$';
    else
        legend_str_vec(ii)="$"+zeta_subtitle+'='+zeta_vec(ii)+'$';
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
        
        xticks=get(ax_mag2,'XTick');
        set(ax_mag2,'XTick',sort(unique([xticks,r_special])));
        set(ax_phase2,'XTick',sort(unique([xticks,r_special])));
    end
end

ylim(ax_phase1,"padded")
ylim(ax_phase2,"padded")
ylim(ax_phase3,"padded")
if N_zeta>1
    figure(Fig_Nyq)
    axis(ax_Lims_Nyq);
    set(gca,'XTick',real_ticks,'YTick',imag_ticks);

    ylim(ax_r,ax_Lims_Nyq(1:2));
    set(ax_r,'YTick',real_ticks);
    ylim(ax_i,ax_Lims_Nyq(3:4));
    set(ax_i,'YTick',imag_ticks);

    ylim(ax_mag1,v(3:4));
end

figure(Fig_3D)
legend(legend_str_vec,'interpreter','latex','Location','bestOutside')

figure(Fig_Nyq)
legend(legend_str_vec,'interpreter','latex','Location','bestOutside')

legend(ax_r,legend_str_vec,'interpreter','latex','Location','northeast')
% legend(ax_i,legend_str_vec,'interpreter','latex','Location','northeast')

legend(ax_mag1,legend_str_vec,'interpreter','latex','Location','northeast')
legend(ax_mag2,legend_str_vec,'interpreter','latex','Location','northeast')
legend(ax_mag3,legend_str_vec,'interpreter','latex','Location','northwest')

set(groot,'DefaultAxesLineStyleOrder','remove')
set(groot,'DefaultLineMarkerSize','remove');