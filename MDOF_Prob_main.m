clearvars
clc
close all

set(groot,'DefaultAxesColorOrder',[0,0,1;0,0,0;1,0,0;0,0.5,0;1,0,1])
set(groot,'DefaultLineLineWidth',1);

f_final=5;
n_points=1000;
n_f_points=1000;
f_column=linspace(0,f_final,n_f_points).';

x_0_col=[-0.05;-0.08];
x_dot_0_col=[-0.01;-0.005];

N=2;
m_vec=[1,1.2];
k_vec=[100,150,100];
% m_vec=12*[1,1,1];
% k_vec=7e3*[0,1,1,1];
% c_vec=10*[1,1,1,1];
% m_vec=[4,5];
% k_vec=[400,500,400];
c_vec=ones(1,N+1);
%c_vec=[8,7,5];
[M_mat,C_mat_temp,K_mat]=N_DOF_sys(m_vec,c_vec,k_vec);
M_mat,K_mat,M_mat\K_mat

m_row=[1,2,2];
n_row=[1,1,2];

%% Undamped system
C_mat=zeros(N,N),t_final=10;t_row=linspace(0,t_final,n_points);MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,true)
filenames=[{'MDOF-Undamped_FRF'};{'MDOF-Undamped_FRF-Nyquist'};{'MDOF-Undamped_IRF'};{'MDOF-Undamped_free_response'};
    {'MDOF-Undamped_Harmonic_response_f1_x1_no_transient'};{'MDOF-Undamped_Harmonic_response_f1_x2_no_transient'};
    {'MDOF-Undamped_Harmonic_response_f1_x1'};{'MDOF-Undamped_Harmonic_response_f1_x2'};
    {'MDOF-Undamped_Harmonic_response_f1_x1-sameY'};{'MDOF-Undamped_Harmonic_response_f1_x2-sameY'};
    {'MDOF-Undamped_Harmonic_response_f2_x1_no_transient'};{'MDOF-Undamped_Harmonic_response_f2_x2_no_transient'};    
    {'MDOF-Undamped_Harmonic_response_f2_x1'};{'MDOF-Undamped_Harmonic_response_f2_x2'};
    {'MDOF-Undamped_Harmonic_response_f2_x1-sameY'};{'MDOF-Undamped_Harmonic_response_f2_x2-sameY'}];
export_figure([2:4],'',filenames(2:4))
export_figure([5:length(filenames)],'||',filenames(5:end))
close all

%To obtain correct phase of undamped FRF
C_mat=10000*eps*ones(N,N),MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,true)
export_figure(1,'',filenames(1))
close all

%% Proportionally viscously damped system
%C_mat=3/4*M_mat+K_mat/300,t_final=6;t_row=linspace(0,t_final,n_points);MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,true)
%C_mat=2/3*M_mat+K_mat/300,t_final=6;t_row=linspace(0,t_final,n_points);MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,true)
C_mat=2/15*M_mat+K_mat/1500,t_final=20;t_row=linspace(0,t_final,n_points);MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,true)
filenames=[{'MDOF-Proportional_FRF'};{'MDOF-Proportional_FRF-Nyquist'};{'MDOF-Proportional_IRF'};{'MDOF-Proportional_free_response'};
    {'MDOF-Proportional_Harmonic_response_f1_x1'};{'MDOF-Proportional_Harmonic_response_f1_x2'};
    {'MDOF-Proportional_Harmonic_response_f1_x1-sameY'};{'MDOF-Proportional_Harmonic_response_f1_x2-sameY'};
    {'MDOF-Proportional_Harmonic_response_f2_x1'};{'MDOF-Proportional_Harmonic_response_f2_x2'};
    {'MDOF-Proportional_Harmonic_response_f2_x1-sameY'};{'MDOF-Proportional_Harmonic_response_f2_x2-sameY'}];
export_figure([1:4],'',filenames(1:4))
export_figure([5:length(filenames)],'||',filenames(5:end))
close all

%% Generally viscously damped system
C_mat=[0.8,-0.1;-0.1,0.3],t_final=15;t_row=linspace(0,t_final,n_points);MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,false)
%C_mat=C_temp,t_final=15;t_row=linspace(0,t_final,n_points);MDOF_Prob(M_mat,C_mat,K_mat,x_0_col,x_dot_0_col,m_row,n_row,t_row,f_column,false)
filenames=[{'MDOF-Viscous_FRF'};{'MDOF-Viscous_FRF-Nyquist'};{'MDOF-Viscous_IRF'};{'MDOF-Viscous_free_response'};
    {'MDOF-Viscous_Harmonic_response_f1_x1'};{'MDOF-Viscous_Harmonic_response_f1_x2'};
    {'MDOF-Viscous_Harmonic_response_f1_x1-sameY'};{'MDOF-Viscous_Harmonic_response_f1_x2-sameY'};
    {'MDOF-Viscous_Harmonic_response_f2_x1'};{'MDOF-Viscous_Harmonic_response_f2_x2'};
    {'MDOF-Viscous_Harmonic_response_f2_x1-sameY'};{'MDOF-Viscous_Harmonic_response_f2_x2-sameY'}];
export_figure([1:4],'',filenames(1:4))
export_figure([5:length(filenames)],'||',filenames(5:end))
close all

%% 
set(groot,'DefaultAxesColorOrder','remove')
set(groot,'DefaultLineLineWidth','remove')