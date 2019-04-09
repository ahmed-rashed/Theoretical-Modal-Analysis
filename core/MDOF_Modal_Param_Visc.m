function [w_r_col, zeta_r_col, w_d_r_col]=MDOF_Modal_Param_Visc(EigValues_vec)

N=length(EigValues_vec)/2;

q_underDamped_vec_temp=find(imag(EigValues_vec)~=0);
q_underDamped_vec=q_underDamped_vec_temp(1:2:end-1);
r_underDamped_vec=(q_underDamped_vec+1)/2;

% q_overDamped_vec_temp=setdiff(1:2*N,q_underDamped_vec_temp);
% q_overDamped_vec_temp=find(imag(EigValues_vec)==0);
q_overDamped_vec=setdiff(1:2:2*N-1,q_underDamped_vec);
% q_overDamped_vec=q_overDamped_vec_temp(1:2:end-1);
r_overDamped_vec=setdiff(1:N,r_underDamped_vec);
% r_overDamped_vec=(q_overDamped_vec+1)/2;

w_r_col=nan(N,1);
zeta_r_col=nan(N,1);
w_d_r_col=nan(N,1);

w_r_col(r_underDamped_vec)=abs(EigValues_vec(q_underDamped_vec));
zeta_r_col(r_underDamped_vec)=-real(EigValues_vec(q_underDamped_vec))./w_r_col(r_underDamped_vec);

zeta_r_col(r_overDamped_vec)=sqrt(1./(1-((EigValues_vec(q_overDamped_vec)-EigValues_vec(q_overDamped_vec+1))./(EigValues_vec(q_overDamped_vec)+EigValues_vec(q_overDamped_vec+1))).^2));
w_r_col(r_overDamped_vec)=-(EigValues_vec(q_overDamped_vec)+EigValues_vec(q_overDamped_vec+1))/2./zeta_r_col(r_overDamped_vec);

w_d_r_col(r_underDamped_vec)=w_r_col(r_underDamped_vec).*sqrt(1-zeta_r_col(r_underDamped_vec).^2);