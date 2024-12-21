function h_m_cols=SDOF_IRF_Visc_mul_m(w_n,zeta_row,t_col)

if ~iscolumn(t_col),error('t_col must be a column vector'),end
if ~isrow(zeta_row),error('zeta_row must be a row vector'),end

N_zeta=size(zeta_row,2);

ind_zeta_lt_1=find(zeta_row<1);
ind_zeta_1=find(zeta_row==1);
ind_zeta_gt_1=setdiff(1:N_zeta,[ind_zeta_lt_1,ind_zeta_1]);

h_m_cols=nan(size(t_col,1),N_zeta);
if ~isempty(ind_zeta_lt_1)
    w_d_row=w_n.*sqrt(1-zeta_row(ind_zeta_lt_1).^2);
    h_m_cols(:,ind_zeta_lt_1)=exp(-zeta_row(ind_zeta_lt_1).*w_n.*t_col).*sin(w_d_row.*t_col)./w_d_row;
end

if ~isempty(ind_zeta_1)
    h_m_cols(:,ind_zeta_1)=exp(-w_n.*t_col).*t_col;
end

if ~isempty(ind_zeta_gt_1)
    beta_row=w_n.*sqrt(zeta_row(ind_zeta_gt_1).^2-1);
    h_m_cols(:,ind_zeta_gt_1)=exp(-zeta_row(ind_zeta_gt_1).*w_n.*t_col).*sinh(beta_row.*t_col)./beta_row;
end

h_m_cols(t_col<0,:)=0;