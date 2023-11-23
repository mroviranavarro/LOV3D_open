function [varargout] = plot_energy_spectra(n_var,m_var,n_v,m_v,energy_s,type,save_name,varargin)
%% Convert into phase and amplitude
energy_spectra_auxR=zeros(size(energy_s));
energy_spectra_amplitude=[];
energy_spectra_phase=[];
spectra_l=size(energy_s,1);
multiple=0; 

for k = 1:length(varargin)
    if strcmpi(varargin{k},'multiple')
        multiple=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end

for i=1:size(energy_s,1)
    k=1; 
    for n_aux=0:max(n_v) %loop over degrees
        for m_aux=0:n_aux %loop over oders
            if m_aux==0
                index=find(m_aux==m_v & n_aux==n_v);
                energy_spectra_auxR(i,index)=energy_s(i,index);
                energy_spectra_amplitude(i,k)=abs(energy_spectra_auxR(i,index));
                energy_spectra_phase(i,k)=0; 
            else
                indexP=find(m_aux==m_v & n_aux==n_v);
                indexN=find(-m_aux==m_v & n_aux==n_v);               
                energy_spectra_auxR(i,indexP)=sqrt(2)*real(energy_s(i,indexP));
                energy_spectra_auxR(i,indexN)=-sqrt(2)*imag(energy_s(i,indexP));                
                energy_spectra_amplitude(i,k)=sqrt(energy_spectra_auxR(i,indexP)^2+energy_spectra_auxR(i,indexN)^2);
                energy_spectra_amplitude(i,k)=abs(energy_s(i,indexP));
                energy_spectra_phase(i,k)=atan2d(-imag(energy_s(i,indexP)),real(energy_s(i,indexP)))/(360);
            end
            n_v2(k)=n_aux;
            m_v2(k)=m_aux;
            k=k+1;
        end
    end
end

%% Make plot
energy_uniform_plot=log10(energy_spectra_amplitude);
energy_uniform_plot(end+1,:)=0;
energy_uniform_plot(:,end+1)=0;
aux=abs(energy_uniform_plot);

fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.45]);
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
tt=tiledlayout(2,1)
ax1=nexttile
h=pcolor(real((energy_uniform_plot)))
set(h, 'EdgeColor', 'white');
colormap(ax1,cmocean('thermal'))
cb=colorbar
if strcmpi(type,'total') 
    cb.Label.String='$\sqrt{{\dot\mathcal{E}_n^m}^2+{\dot\mathcal{E}_n^{-m}}^2}$';
else
    cb.Label.String='$\mid \Delta \dot{e}_n^m \mid/\dot{e}_0^u $';
end
cb.Location='northoutside';
cb.Ticks =[-10:1:0]; 
cb.TickLabels ={'10^{-10}','','10^{-8}','','10^{-6}','','10^{-4}','','10^{-2}','','1',}; 
cb.Label.Interpreter = 'latex';
cb.FontName = 'CMU Serif';
caxis([-10 0]);
i_L=[];
j_l=[];
k=1;
kk=1; 
for i=1:length(n_v2)
        if m_v2(i)==0
           j_l=[j_l k-0.5];
           if n_v2(i)==0
               labels_P{kk}=[''];
           else
                labels_P{kk}=['(' num2str(n_v2(i)) ',' num2str(m_v2(i)) ')'];
           end
          kk=kk+1; 
        elseif m_v2(i)==n_v2(i)
          i_L=[i_L k];
        end
        k=k+1;
end
hold on
plot([2 2],[0 2+spectra_l],'LineWidth',2,'color','w');
hold on
plot([1 1],[0 2+spectra_l],'LineWidth',2,'color','w');
for i=1:length(i_L)
    plot([i_L(i)+1 i_L(i)+1],[0 1+spectra_l],'LineWidth',2,'color','w');
    hold on
end

y_l=[];
for i=1:length(n_var)
    y_l=[y_l i+0.5];
    label_y{i}=['(' num2str(n_var(i)) ',' num2str(m_var(i)) ')'];
    plot([0 length(n_v2)+1],[i i],'LineWidth',2,'color','w');
    hold on 
end

xticks(j_l+1)
xticklabels(labels_P)
yticks(y_l)
yticklabels(label_y)
axis equal
% set-labels y-axis
label_variations_location=[];
label_variations=[];
label_vis_location=[];
ylim([1 1+spectra_l])
set(gca,'fontsize', 30);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
ax1.YAxis.FontSize =16;
ax1.XAxis.FontSize =16;
%xlabel('$(n,m)$','FontSize',25)
%ylabel('$(n_\textrm{LV},m_\textrm{LV})$','FontSize',20)



%% PLOT THE PHASE 
if strcmpi(type,'total') 
else
energy_spectra_phase(end+1,:)=0;
energy_spectra_phase(:,end+1)=0;
ax2=nexttile;
h=pcolor(energy_spectra_phase)
set(h, 'EdgeColor', 'black');
colormap_p_aux=cmocean('curl',1001);
colormap_p_aux=colormap_p_aux(100:end-100,:);
colormap_p=[0 0 0; colormap_p_aux(1:400,:) ; colormap_p_aux(401:end,:) ;0 0 0];
colormap(ax2,colormap_p)
cb2=colorbar;
cb2.Label.String='$\varphi_n^m$';
cb2.Label.Interpreter = 'latex';
cb2.Location='southoutside';
cb2.Label.Interpreter = 'latex';
cb2.FontName = 'CMU Serif';
clim([-0.5 0.5]);

i_L=[];
j_l=[];
k=1;
kk=1; 
for i=1:length(n_v2)
        if m_v2(i)==0
           j_l=[j_l k-0.5];
           if n_v2(i)==0
               labels_P{kk}=[''];
           else
                labels_P{kk}=['(' num2str(n_v2(i)) ',' num2str(m_v2(i)) ')'];
           end
          kk=kk+1; 
        elseif m_v2(i)==n_v2(i)
          i_L=[i_L k];
        end
        k=k+1;
end
hold on
plot([2 2],[0 2+spectra_l],'LineWidth',2,'color','k');
hold on
plot([1 1],[0 2+spectra_l],'LineWidth',2,'color','k');
for i=1:length(i_L)
    plot([i_L(i)+1 i_L(i)+1],[0 1+spectra_l],'LineWidth',2,'color','k');
    hold on
end

y_l=[];
for i=1:length(n_var)
    y_l=[y_l i+0.5];
    label_y{i}=['(' num2str(n_var(i)) ',' num2str(m_var(i)) ')'];
    plot([0 length(n_v2)+1],[i i],'LineWidth',2,'color','k');
    hold on 
end
xticks(j_l+1)
xticklabels('')
yticks(y_l)
yticklabels(label_y)
axis equal
% set-labels y-axis
label_variations_location=[];
label_variations=[];
label_vis_location=[];
ylim([1 1+spectra_l]) 
set(gca,'fontsize', 30);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
ax2.YAxis.FontSize =16;
ax2.XAxis.FontSize =16;
ylabel(tt,'$(n_\textrm{LV},m_\textrm{LV})$','FontSize',25,'interpreter','latex')

tt.TileSpacing = 'tight';
tt.Padding = 'tight';
if isempty(save_name)==0
        export_fig(fig,[save_name],'-png','-opengl','-r300')
end 


end

varargout{2}=energy_spectra_phase; 
varargout{1}=energy_spectra_amplitude;

end