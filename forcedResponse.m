function [y_vec,t_vec]=forcedResponse(h_vec,x_vec, ...
                                        D_t,bRaw)
% If length(h_vec) is less than length(x_vec), h_vec is padded with zeros

if nargin<3
    D_t=1;
else
    if isempty(D_t)
        D_t=1;
    end
end

if nargin<4
    bRaw=false;
else
    if isempty(bRaw)
        bRaw=false;
    end
end

if bRaw
    y_vec=conv(h_vec,x_vec);
    K=length(y_vec);
else
    K=length(x_vec);
    y_vec_temp=conv(h_vec,x_vec);y_vec=y_vec_temp(1:K);
%     y_vec=filter(h_vec,1,x_vec);  %This is equivalent to the above line
end

if sum(x_vec~=0)>1
    y_vec=y_vec*D_t;
end

if nargout>1
    t_vec=(0:K-1)*D_t;
end