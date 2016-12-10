function Hn=MDOF_FRF_Point_Visc(M, C, K, w)

Hn=inv(-w^2*M+i*w*C+K);
