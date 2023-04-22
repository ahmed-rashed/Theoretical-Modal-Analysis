function [w_p_col,zeta_p_col,w_d_p_col]=pole2modal_visc(s_q_pairs_col)
if ~iscolumn(s_q_pairs_col),error('s_q_pairs_col must be a column vector'),end

Q=length(s_q_pairs_col);
P=Q/2;
if P~=round(P),error('s_q_pairs_col must be a vector of even number of elements.'),end

if any(imag(s_q_pairs_col(1:2:end))~=-imag(s_q_pairs_col(2:2:end))),error('Elements of s_q_pairs_col must be consecutive real or complex conjugate pairs.'),end

q_underDamped_vec_temp=find(imag(s_q_pairs_col)~=0);
q_underDamped_vec=q_underDamped_vec_temp(1:2:end-1);
p_underDamped_vec=(q_underDamped_vec+1)/2;

% q_overDamped_vec_temp=setdiff(1:Q,q_underDamped_vec_temp);
% q_overDamped_vec_temp=find(imag(s_q_pairs_col)==0);
q_overDamped_vec=setdiff(1:2:Q-1,q_underDamped_vec);
% q_overDamped_vec=q_overDamped_vec_temp(1:2:end-1);
p_overDamped_vec=setdiff(1:P,p_underDamped_vec);
% p_overDamped_vec=(q_overDamped_vec+1)/2;

w_p_col=nan(P,1);
zeta_p_col=nan(P,1);

w_p_col(p_underDamped_vec)=abs(s_q_pairs_col(q_underDamped_vec));
zeta_p_col(p_underDamped_vec)=-real(s_q_pairs_col(q_underDamped_vec))./w_p_col(p_underDamped_vec);

zeta_p_col(p_overDamped_vec)=sqrt(1./(1-((s_q_pairs_col(q_overDamped_vec)-s_q_pairs_col(q_overDamped_vec+1))./(s_q_pairs_col(q_overDamped_vec)+s_q_pairs_col(q_overDamped_vec+1))).^2));
w_p_col(p_overDamped_vec)=-(s_q_pairs_col(q_overDamped_vec)+s_q_pairs_col(q_overDamped_vec+1))/2./zeta_p_col(p_overDamped_vec);

if nargout>=3
    w_d_p_col=nan(P,1);
    w_d_p_col(p_underDamped_vec)=w_p_col(p_underDamped_vec).*sqrt(1-zeta_p_col(p_underDamped_vec).^2);
end