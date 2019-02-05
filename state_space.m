function state_space()
clc
close all

N=2;

m=[1,1.5];
k=[100,150,100];
c=ones(1,N+1);
[M,CC,K]=N_DOF_sys(m,c,k);

[EigVectors_Normalized, EigValues_vec]=MDOF_Eig_Visc(M, CC, K, false);

[w_r_vec, zeta_r_vec, w_d_r_vec]=MDOF_Modal_Param_Visc(EigValues_vec);

%L1
AA=[zeros(N),-K;-K,-CC];
BB=[-K,zeros(N);zeros(N),M];

A=BB\AA
B=inv(BB)
C=[eye(N),zeros(N)]
D=zeros(N,2*N)

sys=ss(A,B,C,D);
tf(sys)

x_0_col=[-0.05;-0.08];
x_dot_0_col=[-0.01;-0.005];
x_mat_label_col={'x_1';'x_2'};

t_final=30;
n_t=1000;
t_row=linspace(0,t_final,n_t);

h_cols_label_col={'h_{1,1}';'h_{2,1}';'h_{2,2}'};
[h_cols,t_col]=impulse(sys,t_row);
%plot_Forced_Response_Vertically(t_row,h_cols.',h_cols_label_col)

figure
initial(sys,[x_0_col;x_dot_0_col],t_row);

F_0_col=[1;0];
w_F1=[0.5*w_d_r_vec(1),0.9*w_d_r_vec(1),w_d_r_vec(1),1.1*w_d_r_vec(1),(w_d_r_vec(1)+w_d_r_vec(2))/2];
for ii=1:length(w_F1)
    f=F_0_col(1)*sin(w_F1(ii)*t_row.');
    [x_mat,t_col]=lsim(sys,[0*f,0*f,f,0*f],t_row.');
    plot_Forced_Response_Vertically(t_row,x_mat(:,1:N).',x_mat_label_col,f.',{'f_1'},['$w_{0}_1=',num2str(w_F1(ii)/w_d_r_vec(1)),'\ \omega_{\mathrm{d},1}$'])
end

F_0_col=[0;1];
w_F1=[(w_d_r_vec(1)+w_d_r_vec(2))/2,0.9*w_d_r_vec(2),w_d_r_vec(2),1.1*w_d_r_vec(2),2*w_d_r_vec(2)];
for ii=1:length(w_F1)
    f=F_0_col(2)*sin(w_F1(ii)*t_row.');
    [x_mat,t_col]=lsim(sys,[0*f,0*f,0*f,f],t_row.');
    plot_Forced_Response_Vertically(t_row,x_mat(:,1:N).',x_mat_label_col,f.',{'f_2'},['$w_{0}_2=',num2str(w_F1(ii)/w_d_r_vec(2)),'\ \omega_{\mathrm{d},2}$'])
end
