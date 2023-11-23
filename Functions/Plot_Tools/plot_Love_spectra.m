function [varargout] = plot_Love_spectra(n_var,m_var,n_v,m_v,k2t,nf,mf,save_name,varargin)
spectra_l=length(n_var);
spectra_l2=length(n_v);
label_cb=[]; 
x_label=1;
C=-9;
colorbar_p=1;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'label_cb')
        label_cb=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'minimum')
        C=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'colorbar')
        colorbar_p=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'x_label')
        x_label=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end
%% Make plot
k2_log=log10(abs(k2t));
k2_sign=sign(k2t);
k2_plot=k2_sign.*k2_log;
k2_plot=sign(k2t).*log10(1+abs(k2t)/10^C);
k2_plot(end+1,:)=0;
k2_plot(:,end+1)=0;
aux=abs(k2_plot);
fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.4]);
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
h=pcolor(real((k2_plot)))
set(h, 'EdgeColor', 'k');
colormap(cmocean('balance'))

ticks=C:1:0;
k=1;
for j=1:length(ticks)
    for jj=1:1:9
    tick_aux(2*(k-1)+1)=log10(1+jj*10^(ticks(j))/10^C);
    tick_aux(2*(k-1)+2)=-log10(1+jj*10^(ticks(j))/10^C);     
    if jj==1 && mod(ticks(j),2)==0
        TickLabels{2*(k-1)+1}=['$10^{' num2str(ticks(j)) '}$'];
        TickLabels{2*(k-1)+2}=['$-10^{' num2str(ticks(j)) '}$'];
    else
        TickLabels{2*(k-1)+1}=[''];
        TickLabels{2*(k-1)+2}=[''];
    end
    k=k+1;
    end
end
[ticks ticks_ind]=sort(tick_aux);
for i=1:length(ticks_ind)
    TickLabels2{i}=TickLabels{ticks_ind(i)};
end

if colorbar_p==1
cb=colorbar
if isempty(label_cb)
    cb.Label.String=['$\Delta k_{' num2str(nf) ',' num2str(mf) '}^{n,m}/k_{' num2str(nf) '}^u$' ];
else
    cb.Label.String=label_cb;
end

cb.Location='northoutside';
cb.LineWidth= 2;
cb.FontSize=40;
cb.Ticks=ticks; 
cb.Label.Interpreter = 'latex';
set(cb,'TickLabelInterpreter','latex')
cb.Label.FontSize=30;
cb.FontSize=30;
cb.TickLabels =TickLabels2; 
else
    title(['$\Delta k_{' num2str(nf) ',' num2str(mf) '}^{n,m}/k_{' num2str(nf)  '}^u$' ],'FontSize',30);
end
caxis([ticks(1) ticks(end)]);
% leading order
k2tm=k2t;
k2tm(:,1)=0;
%% show leading-order-mode
for i=1:length(n_var)
    [max_k, ind_max]=max(abs(k2tm(i,:)));
    if m_v(ind_max)==0
%         hold on 
%         plot([ind_max ind_max],[i i+1],'LineWidth',1,'color',[0.9 0.9 0.9]);
%         hold on 
%         plot([ind_max+1 ind_max+1],[i i+1],'LineWidth',1,'color',[0.9 0.9 0.9]);
%         hold on
%         plot([ind_max ind_max+1],[i i],'LineWidth',1,'color',[0.9 0.9 0.9]);
%         hold on
%         plot([ind_max ind_max+1],[i+1 i+1],'LineWidth',1,'color',[0.9 0.9 0.9]);
        hold on
        plot([ind_max ind_max+1],[i i+1],'LineWidth',1,'color',[0.9 0.9 0.9]);
        hold on
        plot([ind_max ind_max+1],[i+1 i],'LineWidth',1,'color',[0.9 0.9 0.9]);
        hold on
    else
        k2t_help=abs(k2tm(i,:));
        k2t_help(ind_max)=0;
        [max_k2, ind_max2]=max(abs(k2t_help(:)));
%         hold on 
%         plot([ind_max ind_max],[i i+1],'LineWidth',3,'color','k');
%         hold on 
%         plot([ind_max+1 ind_max+1],[i i+1],'LineWidth',3,'color','k');
%         hold on
%         plot([ind_max ind_max+1],[i i],'LineWidth',3,'color','k');
%         hold on
%         plot([ind_max ind_max+1],[i+1 i+1],'LineWidth',3,'color','k');
%         hold on
          hold on
            plot([ind_max ind_max+1],[i i+1],'LineWidth',1,'color',[0.9 0.9 0.9]);
            hold on
            plot([ind_max ind_max+1],[i+1 i],'LineWidth',1,'color',[0.9 0.9 0.9]);
            hold on
        if abs(max_k2-max_k)<1e-5
%         plot([ind_max2 ind_max2],[i i+1],'LineWidth',3,'color','k');
%         hold on 
%         plot([ind_max2+1 ind_max2+1],[i i+1],'LineWidth',3,'color','k');
%         hold on
%         plot([ind_max2 ind_max2+1],[i i],'LineWidth',3,'color','k');
%         hold on
%         plot([ind_max2 ind_max2+1],[i+1 i+1],'LineWidth',3,'color','k');
%         hold on
        hold on
        plot([ind_max2 ind_max2+1],[i i+1],'LineWidth',1,'color',[0.9 0.9 0.9]);
        hold on
        plot([ind_max2 ind_max2+1],[i+1 i],'LineWidth',1,'color',[0.9 0.9 0.9]);
        hold on
        end
    end
end
    
% set labels x-axis
hold on
i_L=[];
j_l=[];
k=1;
kk=1; 
for i=1:length(n_v)
        if m_v(i)==0
           j_l=[j_l k-0.5];
          labels_P{kk}=['(' num2str(n_v(i)) ',' num2str(m_v(i)) ')'];
          kk=kk+1; 
        elseif m_v(i)==n_v(i)
          i_L=[i_L k];
        end
        k=k+1;
end
plot([2 2],[0 spectra_l2],'LineWidth',2,'color','k');
hold on
plot([1 1],[0 spectra_l2],'LineWidth',2,'color','k');
for i=1:length(i_L)
    plot([i_L(i)+1 i_L(i)+1],[0 spectra_l2],'LineWidth',2,'color','k');
    hold on
end
y_l=[];
hold on
k=1;
for i=1:length(n_var)
    if m_var(i)==0
    y_l=[y_l i+0.5];
    label_y{k}=['(' num2str(n_var(i)) ',' num2str(m_var(i)) ')'];
    plot([0 length(n_v)+1],[i i],'LineWidth',2,'color','k');
    k=k+1;
    end
    hold on 
end
plot([0 length(n_v)+1],[i+1 i+1],'LineWidth',2,'color','k');
xticks(j_l(2:end)+1)
xticklabels(labels_P(2:end))
yticks(y_l)
yticklabels(label_y)
% set-labels y-axis
axis equal
label_variations_location=[];
label_variations=[];
label_vis_location=[];
ylim([1 1+spectra_l]) 


set(gcf,'color','w');  
% set(gca,'fontsize', 22);
 ax=gca;
ax.XAxis.FontSize = 22;
ax.YAxis.FontSize = 22;
set(gcf,'color','w');
ax = gca;
ax.XRuler.Axle.LineWidth = 2;
ax.YRuler.Axle.LineWidth = 2;
% xlabel('Tidal response wavelength ($n,m$)','Fontsize',25)
% ylabel({'Lateral variations ','wavelength ($n$,$m$)'},'Fontsize',25)
if x_label==1
    xlabel('($n,m$)','Fontsize',30)
end
ylabel('$\left(n_\mathrm{LV},m_\mathrm{LV}\right)$','Fontsize',30)
pbaspect([1 1 1])
set(ax,'TickLabelInterpreter','latex')


if isempty(save_name)==0
        export_fig(fig,[save_name '_s'],'-pdf','-opengl','-r200')
end 




end
