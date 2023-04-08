function [w_p_col,zeta_p_col,w_d_p_col]=pole2modal_visc(s_pairs_col)
if ~iscolumn(s_pairs_col),error('s_pairs_col must be a column vector'),end

Q=length(s_pairs_col);
P=Q/2;
if rem(P,1)~=0,error('s_pairs_col must be a vector of even number of elements.'),end

if any(isreal(s_pairs_col)~=repmat(isreal(s_pairs_col(1:2:end)),2)),error('s_pairs_col must consists of consecutive pairs of real or complex numbers.'),end
if any(imag(s_pairs_col(1:2:end))~=-imag(s_pairs_col(2:2:end))),error('Complex numbers in s_pairs_col must be consecutive conjugate pairs.'),end

q_underDamped_vec_temp=find(imag(s_pairs_col)~=0);
q_underDamped_vec=q_underDamped_vec_temp(1:2:end-1);
p_underDamped_vec=(q_underDamped_vec+1)/2;

% q_overDamped_vec_temp=setdiff(1:2*P,q_underDamped_vec_temp);
% q_overDamped_vec_temp=find(imag(EigValues_col)==0);
q_overDamped_vec=setdiff(1:2:2*P-1,q_underDamped_vec);
% q_overDamped_vec=q_overDamped_vec_temp(1:2:end-1);
p_overDamped_vec=setdiff(1:P,p_underDamped_vec);
% p_overDamped_vec=(q_overDamped_vec+1)/2;

w_p_col=nan(P,1);
zeta_p_col=nan(P,1);

w_p_col(p_underDamped_vec)=abs(s_pairs_col(q_underDamped_vec));
zeta_p_col(p_underDamped_vec)=-real(s_pairs_col(q_underDamped_vec))./w_p_col(p_underDamped_vec);

zeta_p_col(p_overDamped_vec)=sqrt(1./(1-((s_pairs_col(q_overDamped_vec)-s_pairs_col(q_overDamped_vec+1))./(s_pairs_col(q_overDamped_vec)+s_pairs_col(q_overDamped_vec+1))).^2));
w_p_col(p_overDamped_vec)=-(s_pairs_col(q_overDamped_vec)+s_pairs_col(q_overDamped_vec+1))/2./zeta_p_col(p_overDamped_vec);

if nargout>=3
    w_d_p_col=nan(P,1);
    w_d_p_col(p_underDamped_vec)=w_p_col(p_underDamped_vec).*sqrt(1-zeta_p_col(p_underDamped_vec).^2);
end