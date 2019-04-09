function MDOF_FRF_ModalSuperposition(f_col,H_cols_SDOF,n_row,m_row, ...
                                                 maxPhaseLag)	 %Optional arguments

if nargin<5
    maxPhaseLag=[];
end

H_cols=sum(H_cols_SDOF,3);
[~,N_FRF,N]=size(H_cols_SDOF);
FigMag=figure;
FigNyq=figure;
for n_FRF=1:N_FRF
    figure(FigMag)
    ax_mag=subplot(4,N_FRF,(n_FRF:N_FRF:n_FRF+2*N_FRF));hold on
    ax_phase=subplot(4,N_FRF,n_FRF+3*N_FRF);hold on

    figure(FigNyq)
    ax_Nyq=subplot(1,N_FRF,n_FRF);hold on

    for r=1:N
        [~,~,h1,h2]=plot_FRF_mag_phase(f_col,H_cols_SDOF(:,n_FRF,r),false,ax_mag,ax_phase,'','',[],maxPhaseLag);

        subplot(ax_Nyq)
        h3=plot_FRF_Nyq(H_cols_SDOF(:,n_FRF,r));
    end
    
    plot_FRF_mag_phase(f_col,H_cols(:,n_FRF),false,ax_mag,ax_phase,'',['H_{',int2str(n_row(n_FRF)),',',int2str(m_row(n_FRF)),'}'],0,maxPhaseLag,'k','LineWidth',1.5*get(h1,'LineWidth'));

    subplot(ax_Nyq)
   %plot_FRF_Nyq(H_cols(:,n_FRF),[],['H_{',int2str(n_row(n_FRF)),',',int2str(m_row(n_FRF)),'}'],0,'k','LineWidth',1.5*get(h3,'LineWidth'));
    plot_FRF_Nyq(H_cols(:,n_FRF),[],['H_{',int2str(n_row(n_FRF)),',',int2str(m_row(n_FRF)),'}'],0,'k');
end