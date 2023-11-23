%% PLOT_ENERGY_MAP
% function used to plot map of energy dissipatoion
%% INPUT
%Energy_Spectra: 
        %Energy_Spectra.n: degrees with non-zero energy
        %Energy_Spectra.m: orders with non-zero energy 
        %Energy_Spectra.n_v: degrees from 0 to Numerics.Nenergy 
        %Energy_Spectra.n_v: orders from 0 to Numerics.Nenergy 
        %Energy_Spectra.energy(radial_point,mode): radial profile of energy spectra
        %Energy_Spectra.energy_integral(mode): radially integrated energy for all non-zero degrees an orders (n,m)
        %Energy_Spectra.energy_integral_v(mode): radially integrated energy for all degrees an orders (n_v,m_v)

% optional variables
    % 'type': 'total' or 'difference'
    % 'save_plot': name of the file where plot is saved 
    % 'limits': specify the limits of the colorscale 
    % 'label': label of the colorbar
    % 'projection': type of projection used for the plot, 'mollweide' is the default and accepts
        % mollweide, 
        % flat 
        % cut: flat plus includes integral along longitude and latitude
        % (used in the paper)
    % 'no_colorbar': can be used to remove colorbar
    % 'no_grid': do not include a grid
    %lon_label: Inlcude longitude label in axis
    %lat_label: Inlcude latitude label in axis
%% OUTPUT 

%%
function [varargout] = plot_energy_map(Energy_Spectra,varargin)
%% TOTAL TIDAL DISSIPATION
n_v=Energy_Spectra.n_v;
m_v=Energy_Spectra.m_v;
energy_s=Energy_Spectra.energy_integral_v;
limits=[];
cut=0; 
projection='mollweide';
label_z=[];
no_colorbar=0; 
no_grid=0; 
type='total';
title_plot='';
save_plot='';
lon_label=1; 
lat_label=1; 
title_color='k';
for k = 1:length(varargin)
    if strcmpi(varargin{k},'limits')
        limits=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'save_plot')
        save_plot=varargin{k+1};  
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'label')
        label_z=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'type')
        type=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'no_colorbar')
        no_colorbar=1; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'no_grid')
        no_grid=1; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'title')
        title_plot=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'title_color')
        title_color=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if  strcmpi(varargin{k},'projection')
        projection=varargin{k+1};
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if  strcmpi(varargin{k},'lon_label')
        lon_label=varargin{k+1};
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if  strcmpi(varargin{k},'lat_label')
        lat_label=varargin{k+1};
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end
%% Get Map
th = linspace(0,pi,500);    % inclination
phi = linspace(-pi,pi,500); % azimuth
[th,phi] = meshgrid(th,phi); 
EnergyV=zeros(size(th)); 
for i=1:length(n_v)
    if abs(energy_s(i))>0
    EnergyV=EnergyV+energy_s(i)*harmonicY(n_v(i),m_v(i),th,phi);
    end
end
aux2=real(EnergyV);
% get colobar limits
if strcmpi(type,'total')
    max_E=max(aux2(:));
    min_E=min(aux2(:)); 
else
    max_E=max(aux2(:));
    min_E=min(aux2(:));  
    max_E2=max([abs(max_E), abs(min_E)]);
end
%% Integrate energy along longitude and latitude 
th_cut=th(1,:); 
phi_cut=phi(:,1); 
Delta_th=th_cut(2)-th_cut(1);
Delta_phi=phi_cut(2)-phi_cut(1);
Energy_Lat=zeros(1,length(th_cut));
Energy_Lon=zeros(1,length(th_cut));
 for i=1:length(th_cut)
     % keep latitude constant, integrate longitude
     Energy_Lat(i)=sum(real(EnergyV(:,i))*Delta_phi);
     % keep longitude constant, integrate latitude
     Energy_Lon(i)=sum(real(EnergyV(i,:)).*sin(th_cut)*Delta_th);
 end
 Energy_Lon=Energy_Lon/2;
 Energy_Lat=Energy_Lat/(2*pi);
varargout{1}=Energy_Lat; 
varargout{2}=Energy_Lon;
sum(Energy_Lat.*sin(th_cut)*Delta_th)/(4*pi);
sum(Energy_Lon*Delta_phi)/(4*pi);
%% Make plot
if strcmpi(projection,'mollweide') 
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.95]);
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(phi*180/pi,90-th*180/pi,real(EnergyV)); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    if strcmpi(type,'total')
        colormap(cmocean('solar',100))
        colormap(cmocean('gray',100))
    else
        colormap(cmocean('balance',100))
    end
    set(0,'defaulttextInterpreter','none')
    if no_grid==0
        m_grid('xaxislocation','middle','xtick',-90:90:180,'ytick',-90:90:90,'fontsize',20,'gridcolor',[0.8 0.8 0.8],'linewidth',3,'linestyle','--','gridcolor','k');
    else
        m_grid('xaxislocation','middle','xtick',-180:30:180,'ytick',-90:15:90,'gridcolor',[0.8 0.8 0.8],'linewidth',3,'linestyle','--','gridcolor','k','fontsize',1);
    end
    if no_colorbar==0
    cb=colorbar('southoutside');
    set(0,'defaulttextInterpreter','latex')
    cb.Label.Interpreter = 'latex';
    if isempty(label_z)==1
        if strcmpi(type,'total')
            cb.Label.String='$\dot{e}$ [-]';
        else
            cb.Label.String='$\Delta\dot{e}/\dot{e}_0^u$ [-]';
        end
    else
        cb.Label.String=label_z; 
        cb.FontSize=50;
    end
    cb.FontName = 'CMU Serif';
    end
    set(gca,'fontsize', 50);
    set(gcf,'color','w');
     if isempty(limits)==0
        caxis([limits(1) limits(2)])
     else
         if strcmpi(type,'total')
            caxis([0 max_E])
        else
           caxis([-max_E2 max_E2])
        end
     end
     if isempty(title_plot)==0
        title(title_plot,'interpreter','latex')
    end
elseif strcmpi(projection,'flat') 
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.95]);
    pcolor(phi*180/pi,90-th*180/pi,real(EnergyV)); 
    shading interp; 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    if strcmpi(type,'total')
        colormap(cmocean('solar',100))
    else
        colormap(cmocean('balance',100))
    end
    set(0,'defaulttextInterpreter','none')
    if no_grid==0
        m_grid('xaxislocation','middle','xtick',[-90:90:180],'ytick',[-90:90:90],'fontsize',20,'gridcolor',[0.8 0.8 0.8],'linewidth',3,'linestyle','--','gridcolor','k');
    else
        m_grid('xaxislocation','middle','xtick',[-90:90:180],'ytick',[-90:90:90],'gridcolor',[0.8 0.8 0.8],'linewidth',3,'linestyle','--','gridcolor','k','fontsize',1);
    end
    if no_colorbar==0
    cb=colorbar('southoutside');
    set(0,'defaulttextInterpreter','latex')
    cb.Label.Interpreter = 'latex';
    if isempty(label_z)==1
        if strcmpi(type,'total')
            cb.Label.String='$\dot{e}$ [-]';
        else
            cb.Label.String='$\Delta\dot{e}$ [-]';
        end
    else
        cb.Label.String=label_z; 
        cb.FontSize=50;
    end
    cb.FontName = 'CMU Serif';
    end
    set(gca,'fontsize', 50);
    set(gcf,'color','w');
     if isempty(limits)==0
        caxis([limits(1) limits(2)])
     else
         if strcmpi(type,'total')
            caxis([0 max_E])
        else
           caxis([-max_E2 max_E2])
        end
     end
     if isempty(title_plot)==0
        title(title_plot,'interpreter','latex')
    end
elseif strcmpi(projection,'cut') 
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.75, 1]);    
    % MAP 
    if no_colorbar==1
        pos1 = [0.1 0.15 0.7 0.53];
    else
        pos1 = [0.1 0.08 0.7 0.6];
    end
    pos2 = [0.8 0.15 0.15 0.53];
    pos3=[0.1 0.68 0.7 0.2];
    subplot('Position',pos1)
    set(gca,'fontsize', 20);
    pcolor(phi*180/pi,90-th*180/pi,real(EnergyV)); 
    shading interp;  
    if strcmpi(type,'total')
        colormap(cmocean('solar',20))
        colormap(cmocean('ice',20))
    else
        colormap(cmocean('balance',21))
    end
    set(0,'defaulttextInterpreter','none')
    if no_colorbar==0
    cb=colorbar('north','Color','w');
    set(0,'defaulttextInterpreter','latex')
    if isempty(label_z)==1
        cb.Label.String='$\Delta\dot{e}$' ;
    else
        cb.Label.String=label_z; 
        cb.FontSize=30;
    end
    cb.Label.Interpreter = 'latex';
    cb.FontName = 'CMU Serif';
    end
    set(gcf,'color','w');
    if isempty(limits)==0
        caxis([limits(1) limits(2)])
     else
         if strcmpi(type,'total')
            caxis([0 max_E])
        else
           caxis([-max_E2 max_E2])
        end
     end    
    xticks(-180:45:180)
    xlim([-180 180]) 
    yticks(-90:45:89)
    ylim([-90 90])
    if lat_label==1
        yticklabels({'$-90^{\circ}$N','$-45^{\circ}$N','$0^{\circ}$','$45^{\circ}$N'})
        ylabel('Latitude','interpreter','latex' )
    else
        yticklabels('')
    end
    if lon_label==1
        xticklabels({'$-180^{\circ}$ (AP)','','$-90^{\circ}$ (LH)','','$0^{\circ}$ (SP)',...
        '','$90^{\circ}$ (TH)','','$0^{\circ}$ (180)'})
        xtickangle(0)
        xlabel('Longitude ','interpreter','latex','FontSize',30)
    else
        xticklabels('')
    end
    set(gca,'fontsize', 30);
    set(gca,'TickLabelInterpreter','latex')
    grid minor
    set(gca, 'layer', 'top');
    text(20,-75,title_plot,'Interpreter','latex','FontSize',40,'color',title_color)
    lon_grid=-90:90:90;
    lat_grid=-45:45:45;
    hold on
    for i=1:length(lon_grid)
        plot([lon_grid(i) lon_grid(i)],[-90 90],'LineWidth',1,'color','k','LineStyle','--')
        hold on
    end
    for i=1:length(lat_grid)
        plot([-180 180],[lat_grid(i) lat_grid(i)],'LineWidth',1,'color','k','LineStyle','--')
        hold on
    end
    
    
    %Longitude averaged
    subplot('Position',pos3)
    box on
    plot(phi_cut*180/pi,Energy_Lon,'LineWidth',3,'color',[0.7 0.7 0.7]);
    xticks(-180:45:180)
    xlim([-180 180]) 
    xticklabels('')
    ax = gca;
    ax.YAxis.FontSize =30;
    ax.XAxis.FontSize =30;
    set(gca, 'XAxisLocation', 'top')
    ylabel('$\frac{1}{{2\dot{e}_0^u}}{\int\Delta\dot{e}\sin\theta d\theta}$','color','k','interpreter','latex','FontSize',20)%,'Units','normalized','Position',[1.05 0.5])
    set(gca, 'YAxisLocation', 'right')
    set(gca,'YTick',[])
    set(gca,'ycolor','k') 
    hold on
    set(gca,'TickLabelInterpreter','latex') 
    min_E=min(Energy_Lon);
    max_E=max(Energy_Lon);
    Delta_E=round((max_E-min_E)/2,2,'significant');
    central_E=round((max_E+min_E)/2,2,'significant');
    y_tick=[central_E-Delta_E central_E central_E+Delta_E];
    ylim([central_E-1.2*Delta_E central_E+1.2*Delta_E  ])
    hold on
    b=max(floor(log10(abs(y_tick))));
    for i=1:length(y_tick)
        a = y_tick(i) * 10^(-b);
        if b==0
            tt=[num2str(a,'%1.2f')];
        else
            tt=[num2str(a,'%1.2f') '$\times 10^{' num2str(b) '}$'];
        end
        text(-178,y_tick(i),tt,'Interpreter','latex','FontSize',30,'VerticalAlignment','middle','HorizontalAlignment','left')
        hold on 
    end 

    yl = ylim;
    for i=1:length(lon_grid)
        plot([lon_grid(i) lon_grid(i)],[yl(1) yl(2)],'LineWidth',1,'color','k','LineStyle','--')
        hold on
    end
    yticks(y_tick)
    set(gcf,'color','w');
    yticklabels('')
    ax.YMinorGrid='off';
    ax.YGrid='on';

    % Latitude averaged 
    subplot('Position',pos2)
    box on 
    plot(Energy_Lat,90-th_cut*180/pi,'LineWidth',3,'color',[0.7 0.7 0.7]);
    set(gca,'TickLabelInterpreter','latex')
    yticks(-90:45:90)
    ylim([-90 90])
    ax = gca;
    ax.YAxis.FontSize =30;
    ax.XAxis.FontSize =30;
    xlabel('$\frac{1}{2\pi\dot{e}_0^u}\int\Delta\dot{e}d\phi$','interpreter','latex','FontSize',20)%,'Units','normalized','Position',[0.63 1.02]);
    ax.XAxisLocation='top';
    yticklabels('')
    hold on
    min_E=min(Energy_Lat);
    max_E=max(Energy_Lat);
    Delta_E=round((max_E-min_E)/2,2,'significant');
    central_E=round((max_E+min_E)/2,2,'significant');
    x_tick=[central_E-Delta_E central_E central_E+Delta_E];
    xlim([central_E-1.4*Delta_E central_E+1.4*Delta_E  ])
    hold on
    b=max(floor(log10(abs(x_tick))));
    for i=1:length(x_tick)
        a = x_tick(i) * 10^(-b);
        text(x_tick(i),82,num2str(a,'%1.2f'),'Interpreter','latex','FontSize',30,'VerticalAlignment','middle','HorizontalAlignment','center')
        hold on 
    end
    hold on 
    if b==0
    else
    text(0.78*x_tick(end),70,['$\times 10^{' num2str(b) '}$'],'Interpreter','latex','FontSize',30,'VerticalAlignment','middle','HorizontalAlignment','center')
    end
    xticklabels('')
    xl = xlim;
    xticks(x_tick)
    for i=1:length(lat_grid)
        plot([xl(1) xl(2)],[lat_grid(i) lat_grid(i)],'LineWidth',1,'color','k','LineStyle','--')
        hold on
    end
    grid minor
    ax.XMinorGrid='off';
    ax.XGrid='on';
 set(gcf,'color','w');
else
    error('Incorrect entry for projection')
end

%% Save Plot
if isempty(save_plot)==0
        %export_fig(fig,save_plot,'-pdf','-opengl')
        export_fig(fig,save_plot,'-transparent','-png','-r300')
        %exportgraphics(fig,[save_plot '.png'], 'Resolution', 300); 
        %f=getframe(gcf);                                                         
        %[X, map] = frame2im(f);                                                                                             
        %imwrite(X,[save_plot '.png']);
        %saveas(fig,[save_plot '.png'])
end
    
end