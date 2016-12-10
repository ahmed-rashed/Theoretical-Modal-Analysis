function [w_r_col, zeta_r_col, w_d_r_col]=MDOF_Modal_Param_Visc(EigValues_mat)

N=size(EigValues_mat,1)/2;

EigValues_vec=diag(EigValues_mat);

UnderDamped_2N_Ind_temp=find(imag(EigValues_vec)~=0);
UnderDamped_2N_Ind=UnderDamped_2N_Ind_temp(1:2:end-1);
UnderDamped_Ind=(UnderDamped_2N_Ind+1)/2;

OverDamped_2N_Ind_temp=find(imag(EigValues_vec)==0);
OverDamped_2N_Ind=OverDamped_2N_Ind_temp(1:2:end-1);
OverDamped_Ind=(OverDamped_2N_Ind+1)/2;

w_r_col=nan(N,1);
zeta_r_col=nan(N,1);
w_d_r_col=nan(N,1);

w_r_col(UnderDamped_Ind)=abs(EigValues_vec(UnderDamped_2N_Ind));
zeta_r_col(UnderDamped_Ind)=-real(EigValues_vec(UnderDamped_2N_Ind))./w_r_col(UnderDamped_Ind);

zeta_r_col(OverDamped_Ind)=sqrt(1./(1-((EigValues_vec(OverDamped_2N_Ind)-EigValues_vec(OverDamped_2N_Ind+1))./(EigValues_vec(OverDamped_2N_Ind)+EigValues_vec(OverDamped_2N_Ind+1))).^2));
w_r_col(OverDamped_Ind)=-(EigValues_vec(OverDamped_2N_Ind)+EigValues_vec(OverDamped_2N_Ind+1))/2./zeta_r_col(OverDamped_Ind);


w_d_r_col(UnderDamped_Ind)=w_r_col(UnderDamped_Ind).*sqrt(1-zeta_r_col(UnderDamped_Ind).^2);