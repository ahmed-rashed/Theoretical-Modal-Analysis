function Hn=MDOF_FRF_Point_Hyst(M, D, K, w)

Hn=inv(K+i*D-w^2*M);
