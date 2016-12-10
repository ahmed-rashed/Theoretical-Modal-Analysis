function MDOF_FRF_main()
clc
close all

set(groot,'DefaultAxesColorOrder',[0,0,1;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultLineMarkerSize',5);
set(groot,'DefaultLineLineWidth',1);

f_final=40;
n_f=8000;
islin=0;
plot_SDOF_FRFs=true;

N=3;

m_vec=10*ones(1,N);
%k_vec=161000/2*ones(1,N+1);
k_vec=[161000/2*ones(1,N),0];
c_vec=8*ones(1,N+1);
%d_vec=2000*ones(1,N+1);

n_row=[1,1,1];
m_row=[1,2,3];

f_column=linspace(0,f_final,n_f).';
w_column=2*pi*f_column;

%Slow FRF calculation
%[M,C,K]=N_DOF_sys(m_vec,c_vec,k_vec);H_w_n_m_cols=MDOF_FRF_slow(@(w)MDOF_FRF_Point_Visc(M, C, K, w), w_column, N, n_row, m_row);
%[M,D,K]=N_DOF_sys(.5*m_vec,d_vec,k_vec);H_w_n_m_cols=MDOF_FRF_slow(@(w)MDOF_FRF_Point_Hyst(M, D, K, w), w_column, N, n_row, m_row);

%Fast FRF calculation
[M,C,K]=N_DOF_sys(m_vec,c_vec,k_vec);
[EigVectors_Normalized, EigValues_mat]=MDOF_Visc_Eig(M, C, K);
H_w_n_m_cols=MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, w_column, n_row, m_row,plot_SDOF_FRFs);

for ii=1:length(n_row)
    H_subtitle=['H_{',int2str(n_row(ii)),',',int2str(m_row(ii)),'}'];
    figure
    plot_FRF_3d(f_column,H_w_n_m_cols(:,ii),'',H_subtitle,1,0);
    
    figure
    ax_r=subplot(2,1,1);
    ax_i=subplot(2,1,2);
    plot_FRF_r_i(f_column,H_w_n_m_cols(:,ii),ax_r,ax_i,'',H_subtitle);
    
    figure
    plot_FRF_Nyq(H_w_n_m_cols(:,ii),H_subtitle);

    figure
    ax_mag=subplot(4,1,[1,2,3]);
    ax_phase=subplot(4,1,4);
    plot_FRF_mag_phase(f_column,H_w_n_m_cols(:,ii),islin,ax_mag,ax_phase,'',H_subtitle,[]);
end

export_figure([1:2],'==',{'MDOF-FRFMag_ModalSuperPos';'FRF-Nyq_ModalSuperPos'})

filenames={'FRF-3D-1';'FRF-RealImag-1';'FRF-Nyq-1';'MDOF-FRFMag1'
           'FRF-3D-2';'FRF-RealImag-2';'FRF-Nyq-2';'MDOF-FRFMag2'
           'FRF-3D-3';'FRF-RealImag-3';'FRF-Nyq-3';'MDOF-FRFMag3'};

export_figure([3:14],'',filenames)

c_vec=8*ones(1,N+1)*10;
%Fast FRF calculation
[M,C,K]=N_DOF_sys(m_vec,c_vec,k_vec);
[EigVectors_Normalized, EigValues_mat]=MDOF_Visc_Eig(M, C, K);
MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, w_column, n_row, m_row,plot_SDOF_FRFs);
export_figure([15:16],'==',{'MDOF-FRFMag_ModalSuperPos1';'FRF-Nyq_ModalSuperPos1'})

set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineMarkerSize','remove');
set(groot,'DefaultLineLineWidth','remove');