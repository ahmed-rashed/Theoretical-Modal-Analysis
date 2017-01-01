function H_w_n_m_cols= ...
MDOF_FRF_Visc(EigValues_mat, EigVectors_Normalized, w_column, n_row, m_row, ...
		      plot_SDOF_FRFs,maxPhaseLag)	 %Optional arguments

if nargin<6
    plot_SDOF_FRFs=false;
else
    if isempty(plot_SDOF_FRFs)
        plot_SDOF_FRFs=false;
    end
end

if nargin<7
    maxPhaseLag=[];
end

N=size(EigVectors_Normalized,1);
n_f=size(w_column,1);
if (any(size(n_row)~=size(m_row)));error('Dimensions of n_row and m_row must be identical');end
i_col=size(n_row,2);

if plot_SDOF_FRFs
    ax_mag=zeros(1,i_col);
    ax_phase=zeros(1,i_col);
    ax_Nyq=zeros(1,i_col);
    FigMag=figure;
    FigNyq=figure;
    for ii=1:i_col
        figure(FigMag)
        ax_mag(ii)=subplot(4,i_col,[ii:i_col:ii+2*i_col]);hold on
        ax_phase(ii)=subplot(4,i_col,ii+3*i_col);hold on
        
        figure(FigNyq)
        ax_Nyq(ii)=subplot(1,i_col,ii);hold on
    end
end

H_w_n_m_cols_SDOF=zeros(n_f,i_col);
H_w_n_m_cols=zeros(n_f,i_col);
A_ind_row=sub2ind([N,N],n_row,m_row);
for r=1:2*N
    A_r=EigVectors_Normalized(:,r)*EigVectors_Normalized(:,r).';
    A_r_temp_row=A_r(A_ind_row);
    H_w_n_m_cols_SDOF=H_w_n_m_cols_SDOF+(1./(1i*w_column-EigValues_mat(r,r)))*A_r_temp_row;
    
    if imag(EigValues_mat(r,r))~=0 && mod(r,2)~=0   %complex eigenvalue and odd r
        continue
    else     %real eigenvalue or even r
        if plot_SDOF_FRFs
            for ii=1:i_col
                [ax_mag(ii),ax_phase(ii),h1,h2]=plot_FRF_mag_phase(w_column/2/pi,H_w_n_m_cols_SDOF(:,ii),false,ax_mag(ii),ax_phase(ii),'','',[],maxPhaseLag);
                
                figure(FigNyq)
                subplot(ax_Nyq(ii))
                h3=plot_FRF_Nyq(H_w_n_m_cols_SDOF(:,ii));
            end
        end
    end
    
    H_w_n_m_cols=H_w_n_m_cols+H_w_n_m_cols_SDOF;
    H_w_n_m_cols_SDOF=zeros(n_f,i_col);
end

if plot_SDOF_FRFs
    for ii=1:i_col
        plot_FRF_mag_phase(w_column/2/pi,H_w_n_m_cols(:,ii),false,ax_mag(ii),ax_phase(ii),'',['H_{',int2str(n_row(ii)),',',int2str(m_row(ii)),'}'],0,maxPhaseLag,'k','LineWidth',1.5*get(h1,'LineWidth'));

        subplot(ax_Nyq(ii))
        %plot_FRF_Nyq(H_w_n_m_cols(:,ii),['H_{',int2str(n_row(ii)),',',int2str(m_row(ii)),'}'],0,'k','LineWidth',1.5*get(h3,'LineWidth'));
        plot_FRF_Nyq(H_w_n_m_cols(:,ii),['H_{',int2str(n_row(ii)),',',int2str(m_row(ii)),'}'],0,'k');
    end
end