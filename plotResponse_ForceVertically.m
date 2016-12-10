function plotResponse_ForceVertically(t_row,x_rows,x_rows_label_col,  ...
                        f_rows,f_rows_labels_col,axisTitle_latex,sameScale_y1) %Optional arguments

if nargin<4
    f_rows=[];
end

if nargin<7
    sameScale_y1=false;
end
if sameScale_y1
    axisTitle_latex=[axisTitle_latex,'; (same ',x_rows_label_col{1},' limits)'];
end

N_signals=size(x_rows,1);
%figure
yLimitsMin=inf;
yLimitsMax=-inf;
AX1=zeros(N_signals,1);
for n=1:N_signals
    subplot(N_signals,1,n)
    if nargin<4 || isempty(f_rows)
        plot(t_row,x_rows(n,:));
        AX1(n)=gca;
        
        if n~=N_signals
            set(AX1(n),'XTickLabel',[]);
        end
    else
        if size(f_rows,1)==1
            [AX,h1,h2]=plotyy(t_row,x_rows(n,:),t_row,f_rows);
        else
            [AX,h1,h2]=plotyy(t_row,x_rows(n,:),t_row,f_rows(n,:));
        end
        set(h1,'Color','r')
        set(AX(1),'YColor','r')
        set(h2,'Color','b')
        set(AX(2),'YColor','b')
        uistack(AX(1),'top')
        colorTemp=AX(1).Color;
        boxTemp=AX(1).Box;
        visibleTemp=AX(1).XAxis.Visible;
        AX(1).Color=AX(2).Color;
        AX(2).Color=colorTemp;
        AX(1).Box='off';
        AX(2).Box='off';
        AX(1).XAxis.Visible='on';
        AX(2).XAxis.Visible='on';
        AX(2).XAxisLocation='top'
        set(AX(2),'XTickLabel',[]);
        
        AX1(n)=AX(1);
        if length(f_rows_labels_col)==1
            if ~isempty(f_rows_labels_col{1})
                ylabel(f_rows_labels_col{1},'interpreter','latex');
            end
        elseif ~isempty(f_rows_labels_col)
            if ~isempty(f_rows_labels_col(n))
                ylabel(AX(2),f_rows_labels_col{n},'interpreter','latex');
            end
        end
        set(h1,'LineWidth',1.5*get(h2,'LineWidth'));
        
        %axis(AX(2),'tight');
        ylim(AX(2),max(abs(ylim(AX(2))))*[-1,1]);
        %set(AX(2),'YTickMode','auto');

        if n~=N_signals
            for ss=1:2
                set(AX(ss),'XTickLabel',[]);
            end
        end
    end
    
    if length(x_rows_label_col)==1
        if ~isempty(x_rows_label_col{1})
            ylabel(x_rows_label_col{1},'interpreter','latex')
        end
    elseif ~isempty(x_rows_label_col)
        if ~isempty(x_rows_label_col{n})
            ylabel(x_rows_label_col{n},'interpreter','latex')
        end
    end
    
    axis(AX1(n),'tight');

    ylim(AX1(n),max(abs(ylim(AX1(n))))*[-1,1]);
    set(AX1(n),'YTickMode','auto');

    if sameScale_y1 && ~all(isnan(x_rows(n,:)))
        yLimits=get(AX1(n),'ylim');
        yLimitsMin=min(yLimitsMin,yLimits(1));
        yLimitsMax=max(yLimitsMax,yLimits(2));
    end
    
    %xlabel('$t$','interpreter','latex')
end
xlabel('$t$','interpreter','latex')

if sameScale_y1 && yLimitsMin~=inf && yLimitsMax~=-inf
    for n=1:N_signals
        ylim(AX1(n),[yLimitsMin,yLimitsMax]);
        set(AX1(n),'YTickMode','auto')
    end
end

if nargin>5 && ~isempty(axisTitle_latex)
    title(AX1(1),axisTitle_latex,'interpreter','latex')
end