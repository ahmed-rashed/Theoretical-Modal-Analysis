function plot_Forced_Response_Vertically(t_row,x_rows, ...
                        f_rows,title_str,f_title_str,r_str_col,x_title_str,sameScale_y1,sameScale_y2) %Optional arguments

if nargin<3
    f_rows=[];
end

if nargin<8
    sameScale_y1=true;
end

if nargin<9
    sameScale_y2=true;
end

if ~sameScale_y1
    x_title_str=x_title_str+'; (different limits)';
end

if ~sameScale_y2
    f_title_str=f_title_str+'; (different limits)';
end

N_signals=size(x_rows,1);
%figure
oAx_vec=gobjects(N_signals,1);
tile1=tiledlayout(N_signals,1,"TileSpacing","tight");
for n=1:N_signals
    oAx_vec(n)=nexttile;
    if nargin<3 || isempty(f_rows)
        plot(t_row,x_rows(n,:));
        
        if n~=N_signals
            oAx_vec(n).XTickLabel=[];
        end
    else
        yyaxis right
        if size(f_rows,1)==1
            h2=plot(t_row,f_rows);
        else
            h2=plot(t_row,f_rows(n,:));
        end
        if isscalar(r_str_col)
            if r_str_col(1)~=""
                ylabel(r_str_col(1),'interpreter','latex');
            end
        elseif ~isempty(r_str_col)
            if r_str_col(n)~=""
                ylabel(r_str_col(n),'interpreter','latex');
            end
        end
        ylim(max(abs(ylim))*[-1,1]);

        yyaxis left
        if all(isnan(x_rows(n,:)))
            ylim(100*eps*[-1,1]);
        else
            h1=plot(t_row,x_rows(n,:));
            set(h1,'LineWidth',3*get(h2,'LineWidth'));
        end

        if n~=N_signals
            set(oAx_vec(n),'XTickLabel',[]);
        end
    end
    
    if ~isempty(get(oAx_vec(n),'Children'))
        axis('tight');
        ylim(max(abs(ylim))*[-1,1]);
    end
end
xlabel('$t$','interpreter','latex')

if sameScale_y1
    linkaxes(oAx_vec,'y');
end

if sameScale_y2
    for n=1:N_signals
        yyaxis(oAx_vec(n),"right");
    end
    linkaxes(oAx_vec,'y');
end

if nargin>4 && ~isempty(title_str)
    title(tile1,title_str,'interpreter','latex')
end

ax_hidden=axes(tile1,'Visible','off');
ax_hidden.Layout.TileSpan=[N_signals,1];

yyaxis(ax_hidden,'left');
ylabel([x_title_str;""],'interpreter','latex','FontSize',12,'Visible','on');

yyaxis(ax_hidden,'right');
ylabel(["";"";f_title_str],'interpreter','latex','FontSize',12,'Visible','on');