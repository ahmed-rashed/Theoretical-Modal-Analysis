function h_cols= ...
MDOF_IRF_Visc(EigValues_vec,EigVectors_Normalized,t_col,m_vec,n_vec)

[N,Q]=size(EigVectors_Normalized);
if Q~=2*N,error('EigVectors_Normalized should be N x 2N matrix.'),end
K=size(t_col,1);
if ~isvector(m_vec),error('m_vec must be a vector!'),end
if (any(size(m_vec)~=size(n_vec)));error('Dimensions of m_row and n_row must be identical');end
N_IRF=length(m_vec);

h_cols=zeros(K,N_IRF);
A_ind_row=sub2ind([N,N],m_vec(:).',n_vec(:).');

for r=1:2*N
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    A_r_temp_row=A_r(A_ind_row);
    h_cols=h_cols+exp(EigValues_vec(r)*t_col)*A_r_temp_row;
    
% %For Display only
% if imag(EigValues_vec(r))~=0 && mod(r,2)~=0
%     disp('IRF Response parameters');
%     r
%     2*abs(A_r)
% end
end

ind=find(abs(imag(h_cols))<=10*eps);
h_cols(ind)=real(h_cols(ind));