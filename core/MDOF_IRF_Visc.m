function h_cols= ...
MDOF_IRF_Visc(EigValues_vec,EigVectors_Normalized,t_col,m_vec,n_vec)

[P,Q]=size(EigVectors_Normalized);
if Q~=2*P,error('EigVectors_Normalized should be P x 2N matrix.'),end
K=size(t_col,1);
if ~isvector(m_vec),error('m_vec must be a vector!'),end
if (any(size(m_vec)~=size(n_vec)));error('Dimensions of m_row and n_row must be identical');end
N_IRF=length(m_vec);

h_cols=zeros(K,N_IRF);
A_ind_row=sub2ind([P,P],m_vec(:).',n_vec(:).');

for q=1:Q
    A_q=EigVectors_Normalized(:,q)*EigVectors_Normalized(:,q).';
    A_q_temp_row=A_q(A_ind_row);
    h_cols=h_cols+exp(EigValues_vec(q)*t_col)*A_q_temp_row;
    
% %For Display only
% if imag(EigValues_vec(q))~=0 && mod(q,2)~=0
%     disp('IRF Response parameters');
%     q
%     2*abs(A_q)
% end
end

ind=find(abs(imag(h_cols))<=10*eps);
h_cols(ind)=real(h_cols(ind));