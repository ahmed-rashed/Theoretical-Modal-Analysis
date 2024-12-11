function [s_q_pairs_col,A_q_mat_pages]=modal2pole_visc(w_p_col,zeta_p_col,A_p_mat_pages)

if ~iscolumn(w_p_col),error('w_p_col must be a column vector'),end
P=length(w_p_col);
if any(size(zeta_p_col)~=[P,1]),error('w_p_col and zeta_p_col must have identical sizes!'),end

Q=2*P;

p_underDamped_vec=find(zeta_p_col<1);
p_overDamped_vec=setdiff(1:P,p_underDamped_vec);

s_q_pairs_col=nan(Q,1);
w_d_underDamped_col=w_p_col(p_underDamped_vec).*sqrt(1-zeta_p_col(p_underDamped_vec).^2);
s_q_pairs_col(2*p_underDamped_vec-1)=-zeta_p_col(p_underDamped_vec).*w_p_col(p_underDamped_vec)+1i*w_d_underDamped_col;
s_q_pairs_col(2*p_underDamped_vec)=conj(s_q_pairs_col(2*p_underDamped_vec-1));

s_q_pairs_col(2*p_overDamped_vec-1)=-zeta_p_col(p_overDamped_vec).*w_p_col(p_overDamped_vec)+w_p_col(p_overDamped_vec).*sqrt(zeta_p_col(p_overDamped_vec).^2-1);
s_q_pairs_col(2*p_overDamped_vec)=-zeta_p_col(p_overDamped_vec).*w_p_col(p_overDamped_vec)-w_p_col(p_overDamped_vec).*sqrt(zeta_p_col(p_overDamped_vec).^2-1);

if (nargin>2) && (nargout>1)
    A_q_mat_pages=repmat(A_p_mat_pages,1,1,2);
    A_q_mat_pages(:,:,2*p_underDamped_vec)=conj(A_p_mat_pages(:,:,p_underDamped_vec));
    A_q_mat_pages(:,:,2*p_overDamped_vec)=-A_p_mat_pages(:,:,p_overDamped_vec);
end