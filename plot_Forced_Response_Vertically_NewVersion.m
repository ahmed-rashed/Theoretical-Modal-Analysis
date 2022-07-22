function plot_Forced_Response_Vertically(t_row,x_rows,x_rows_Latex_sym_col, ...
                        f_rows,f_rows_label_col,axisTitle_latex,sameScale_y1) %Optional arguments

if nargin<4
    f_rows=[];
end

if nargin<7
    sameScale_y1=false;
end
if sameScale_y1
    axisTitle_latex=axisTitle_latex+'; (same $'+x_rows_Latex_sym_col+'$ limits)';
end

N_signals=size(x_rows,1);
%figure
yLimitsMin=inf;
yLimitsMax=-inf;
for n=1:N_signals
    AX=subplot(N_signals,1,n);
    if n~=N_signals
         AX.XTickLabel=[];
    end

    if nargin<4 || isempty(f_rows)
        plot(t_row,x_rows(n,:));
    else
        if size(f_rows,1)==1
            yyaxis right;
            plot(t_row,f_rows,'b');
            AX.YColor='b';
            
            yyaxis left;
            plot(t_row,f_rows,'r','LineWidth',1.5*get(groot,'DefaultLineLineWidth'));
            AX.YColor='r';
        else
            yyaxis right;
            plot(t_row,f_rows(n,:),'b');
            AX.YColor='b';
            
            yyaxis left;
            plot(t_row,x_rows(n,:),'r','LineWidth',1.5*get(groot,'DefaultLineLineWidth'));
            AX.YColor='r';
        end
        uistack(AX.YAxis(1),'top')
        colorTemp=AX.YAxis(1).Color;
        AX.YAxis(1).Color=AX.YAxis(2).Color;
        AX.YAxis(2).Color=colorTemp;

        
        yyaxis right;
        if length(f_rows_label_col)==1
            if ~isempty(f_rows_label_col{1})
                ylabel(f_rows_label_col{1},'interpreter','latex');
            end
        elseif ~isempty(f_rows_label_col)
            if ~isempty(f_rows_label_col(n))
                ylabel(f_rows_label_col{n},'interpreter','latex');
            end
        end
    end
    
    yyaxis left;
    if length(x_rows_Latex_sym_col)==1
        if ~isempty(x_rows_Latex_sym_col{1})
            ylabel(x_rows_Latex_sym_col{1},'interpreter','latex')
        end
    elseif ~isempty(x_rows_Latex_sym_col)
        if ~isempty(x_rows_Latex_sym_col{n})
            ylabel(x_rows_Latex_sym_col{n},'interpreter','latex')
        end
    end
    
    axis tight

    yyaxis left;
    ylim(max(abs(ylim))*[-1,1]);
    AX.YTickMode='auto';

    if nargin<4 || isempty(f_rows)
        yyaxis right;
        ylim(max(abs(ylim))*[-1,1]);
        AX.YTickMode='auto';
    end

    if sameScale_y1 && ~all(isnan(x_rows(n,:)))
        yyaxis left;
        yLimits=ylim;
        yLimitsMin=min(yLimitsMin,yLimits(1));
        yLimitsMax=max(yLimitsMax,yLimits(2));
    end
    
    %xlabel('$t$','interpreter','latex')
end
xlabel('$t$','interpreter','latex')

if sameScale_y1 && yLimitsMin~=inf && yLimitsMax~=-inf
    for n=1:N_signals
        AX=subplot(N_signals,1,n);
        yyaxis left;
        ylim([yLimitsMin,yLimitsMax]);
        AX.YTickMode='auto';
    end
end

if nargin>5 && ~isempty(axisTitle_latex)
    subplot(N_signals,1,1)
    title(axisTitle_latex,'interpreter','latex')
end