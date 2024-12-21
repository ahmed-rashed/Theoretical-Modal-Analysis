function y_step_cols=SDOF_Vehicle_Step_Response(Y_0,w_n,zeta_row,t_col)

if ~iscolumn(t_col),error('t_col must be a column vector'),end
if ~isrow(zeta_row),error('zeta_row must be a row vector'),end

N_zeta=size(zeta_row,2);

ind_zeta_lt_1=find(zeta_row<1);
ind_zeta_1=find(zeta_row==1);
ind_zeta_gt_1=setdiff(1:N_zeta,[ind_zeta_lt_1,ind_zeta_1]);

y_step_cols=nan(size(t_col,1),N_zeta);
if ~isempty(ind_zeta_lt_1)
    w_d_row=w_n.*sqrt(1-zeta_row(ind_zeta_lt_1).^2);
    y_step_cols(:,ind_zeta_lt_1)=1+exp(-zeta_row(ind_zeta_lt_1).*w_n.*t_col).*(zeta_row(ind_zeta_lt_1).*w_n.*sin(w_d_row.*t_col)/w_d_row-cos(w_d_row.*t_col));
end

if ~isempty(ind_zeta_1)
    y_step_cols(:,ind_zeta_1)=1+exp(-w_n.*t_col).*(w_n.*t_col-1);
end
if ~isempty(ind_zeta_gt_1)
    beta_row=w_n.*sqrt(zeta_row(ind_zeta_gt_1).^2-1);
    y_step_cols(:,ind_zeta_gt_1)=1+exp(-zeta_row(ind_zeta_gt_1).*w_n.*t_col).*(zeta_row(ind_zeta_gt_1).*w_n.*sinh(beta_row.*t_col)/beta_row-cosh(beta_row.*t_col));
end

y_step_cols(t_col<0,:)=0;
y_step_cols=Y_0.*y_step_cols;