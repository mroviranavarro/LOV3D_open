%% GET_STRESS_STRAIN_MAP
% AUTHOR: M. Rovira-Navarro 
% USE: plot response in map
%% INPUT
    % y_LatLon: solution in map
        % y_LatLon.nf: degree of the forcing
        % y_LatLon.mf: order of the forcing 
        % y_LatLon.lon: longitude
        % y_LatLon.lat: latitude 
        % y_LatLon.r: radial point
        % y.y(longitude,latitude,radial_point,X) 
            %y.y(longitude,latitude,radial_point,1):  Gravitational Potential 
            %y.y(longitude,latitude,radial_point,2):  Displacement e_r component 
            %y.y(longitude,latitude,radial_point,3):  Displacement e_theta component 
            %y.y(longitude,latitude,radial_point,4):  Displacement e_phi component
            %y.y(longitude,latitude,radial_point,5):  stress   e_r e_r component 
            %y.y(longitude,latitude,radial_point,6):  stress   e_r e_theta component 
            %y.y(longitude,latitude,radial_point,7):  stress   e_r e_phi component 
            %y.y(longitude,latitude,radial_point,8):  stress   e_theta e_r compone
            %y.y(longitude,latitude,radial_point,9):  stress   e_theta e_theta component 
            %y.y(longitude,latitude,radial_point,10): stress   e_theta e_phi component 
            %y.y(longitude,latitude,radial_point,11): stress   e_phi e_r component 
            %y.y(longitude,latitude,radial_point,12): stress   e_phi e_theta component 
            %y.y(longitude,latitude,radial_point,13): stress   e_phi e_phi component 
            %y.y(longitude,latitude,radial_point,14): strain   e_r e_r component 
            %y.y(longitude,latitude,radial_point,15): strain   e_r e_theta component 
            %y.y(longitude,latitude,radial_point,16): strain   e_r e_phi component 
            %y.y(longitude,latitude,radial_point,17): strain   e_theta e_r compone
            %y.y(longitude,latitude,radial_point,18): strain   e_theta e_theta component 
            %y.y(longitude,latitude,radial_point,19): strain   e_theta e_phi component 
            %y.y(longitude,latitude,radial_point,20): strain   e_phi e_r component 
            %y.y(longitude,latitude,radial_point,21): strain   e_phi e_theta component 
            %y.y(longitude,latitude,radial_point,22): strain   e_phi e_phi component 
            % optional variables 
        % save_plot: save plot to file
        % radial_point: specify the radial poin 
        % field_name: field to be plotted
            % all (default) 
            % grav_potential
            % displacements
            % stress (all components)
            % strain (all components)
            % stress_rr
            % stress_rtheta
            % stress_rphi
            % stress_thetar
            % stress_thetaphi
            % stress_phir
            % stress_phitheta
            % stress_phiphi
            % strain_rr
            % strain_rtheta
            % strain_rphi
            % strain_thetar
            % strain_thetaphi
            % strain_phir
            % strain_phitheta
            % strain_phiphi
            % test_strain_stress           
%% OUTPUT 
    % map 
function [] = plot_map(y_LatLon,Interior_Model,varargin)
%% OPTIONAL INPUTS
save_plot=0; 
field_name='all';
r_point=length(y_LatLon.r)-1; 
plot_title='';
for k = 1:length(varargin)
    if strcmpi(varargin{k},'field_name')
        field_name=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'radial_point')
        r_point=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
     if strcmpi(varargin{k},'save_plot')
        save_name=varargin{k+1};
        save_plot=1; 
        varargin{k+1}=[]; 
        varargin{k}=[];
     end
    if strcmpi(varargin{k},'plot_title')
        plot_title=varargin{k+1};
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end
%% MAKE PLOT
switch field_name
    case 'all' %%%------------ ALL COMPONENTS 
        fig=figure;
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.9, 0.9]);
        
        p = uipanel('Parent',fig,'BorderType','none'); 
        if isempty(plot_title)==0
        p.Title = plot_title; 
        p.TitlePosition = 'centertop'; 
        p.FontSize = 18;
        p.FontWeight = 'bold';
        end

        tit={'$\Re(\phi)$','$\Re(u)$','$\Re(\nabla\cdot u)$','$\Re(\phi^F)$','$\Re(\mu)$','$Im(\mu)$','$\Re(\sigma_{rr})$','$\Re(\sigma_{r\theta})$','$\Re(\sigma_{r\phi})$','$\Re(\sigma_{\theta r})$','$\Re(\sigma_{\theta \theta})$','$\Re(\sigma_{\theta \phi})$','$\Re(\sigma_{\phi r})$','$\Re(\sigma_{\phi \theta})$','$\Re(\sigma_{\phi \phi})$',...
            '$\Re(\epsilon_{rr})$','$\Re(\epsilon_{r\theta})$','$\Re(\epsilon_{r\phi})$','$\Re(\epsilon_{\theta r})$','$\Re(\epsilon_{\theta \theta})$','$\Re(\epsilon_{\theta \phi})$','$\Re(\epsilon_{\phi r})$','$\Re(\epsilon_{\phi \theta})$','$\Re(\epsilon_{\phi \phi})$'};
        plot_n=[1,2,3,4,5,6,7,8,9,13,14,15,19,20,21,10,11,12,16,17,18,22,23,24];
        plot_c=[1,2,3,4,4,4,5:1:22];
        for j=1:length(plot_n)
        subplot(4,6,plot_n(j),'Parent',p)
        if j==3
            to_plot=real(y_LatLon.y(:,:,r_point,14)+y_LatLon.y(:,:,r_point,18)+y_LatLon.y(:,:,r_point,22));
        elseif j==4
            to_plot=real(y_LatLon.forcing);
        elseif j==5
            to_plot=real(y_LatLon.mu);
        elseif j==6
            to_plot=imag(y_LatLon.mu);
        else
            to_plot=real(y_LatLon.y(:,:,r_point,plot_c(j)));
        end
        h=pcolor(y_LatLon.lon*90/pi,y_LatLon.lat*90/pi,to_plot);
        set(h, 'EdgeColor', 'none');
        if j==2
            dis_theta=(real(y_LatLon.y(1:10:end,1:10:end,r_point,3)));
            dis_phi=(real(y_LatLon.y(1:10:end,1:10:end,r_point,4)));
            hold on
            quiver(y_LatLon.lon(1:10:end,1:10:end)*90/pi,y_LatLon.lat(1:10:end,1:10:end)*90/pi,dis_phi,-dis_theta,'k','LineWidth',2)
        end
        colorbar
        if j==5  || j==6 
            colormap(cmocean('balance',20))
            c_lim=max(abs(to_plot(:))-1);
            if c_lim>0
            caxis([1-c_lim 1+c_lim])
            end
        else
            colormap(cmocean('balance',20))
            c_lim=max(abs(to_plot(:)));
            if c_lim>0
            caxis([-c_lim c_lim])
            else
            end
        end
        title(tit{j},'interpreter','latex')
        set(gca,'Box','on');
        set(gca,'fontsize', 18);
        end
    case 'grav_potential' %%%------------ Gravitational potential 
        fig=figure;
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.5]);
        to_plot=real(y_LatLon.y(:,:,r_point,1));
        h=pcolor(y_LatLon.lon*90/pi,y_LatLon.lat*90/pi,to_plot);
        set(h, 'EdgeColor', 'none');
        colorbar
        colormap(cmocean('balance',20))
        c_lim=max(abs(to_plot(:)));
        caxis([-c_lim c_lim])
        title('$\Re(\phi)$','interpreter','latex')
        set(gca,'fontsize', 18);
    case 'displacements' %%%------------ Displacements 
        fig=figure;
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.5]);
        to_plot=real(y_LatLon.y(:,:,r_point,2));
        h=pcolor(y_LatLon.lon*90/pi,y_LatLon.lat*90/pi,to_plot);
        set(h, 'EdgeColor', 'none');
        dis_theta=(real(y_LatLon.y(1:10:end,1:10:end,r_point,3)));
        dis_phi=(real(y_LatLon.y(1:10:end,1:10:end,r_point,4)));
        hold on
        quiver(y_LatLon.lon(1:10:end,1:10:end)*90/pi,y_LatLon.lat(1:10:end,1:10:end)*90/pi,dis_phi,-dis_theta,'k','LineWidth',2)
        colorbar
        colormap(cmocean('balance',20))
        c_lim=max(abs(to_plot(:)));
        caxis([-c_lim c_lim])
        title('$u$','interpreter','latex')
        set(gca,'Box','on');
        set(gca,'fontsize', 18);    
    case 'stress' %%%------------ Stress 
        fig=figure;
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.9]);
        tit={'$\Re(\sigma_{rr})$','$\Re(\sigma_{r\theta})$','$\Re(\sigma_{r\phi})$','$\Re(\sigma_{\theta r})$','$\Re(\sigma_{\theta \theta})$','$\Re(\sigma_{\theta \phi})$','$\Re(\sigma_{\phi r})$','$\Re(\sigma_{\phi \theta})$','$\Re(\sigma_{\phi \phi})$'};
        plot_c=[5:1:13];
        for i=1:9
            to_plot=real(y_LatLon.y(:,:,r_point,plot_c(i)));
            subplot(3,3,i)
            h=pcolor(y_LatLon.lon*90/pi,y_LatLon.lat*90/pi,to_plot);
            set(h, 'EdgeColor', 'none');
            colorbar
            colormap(cmocean('balance',20))
            c_lim=max(abs(to_plot(:)));
            caxis([-c_lim c_lim])
            title(tit{i},'interpreter','latex')
            set(gca,'Box','on');
            set(gca,'fontsize', 18); 
        end
    case 'strain'  %%%------------ Strain 
        fig=figure;
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.9]);
        tit={'$\Re(\epsilon_{rr})$','$\Re(\epsilon_{r\theta})$','$\Re(\epsilon_{r\phi})$','$\Re(\epsilon_{\theta r})$','$\Re(\epsilon_{\theta \theta})$','$\Re(\epsilon_{\theta \phi})$','$\Re(\epsilon_{\phi r})$','$\Re(\epsilon_{\phi \theta})$','$\Re(\epsilon_{\phi \phi})$'};
        plot_c=[14:1:22];
        for i=1:9
            to_plot=real(y_LatLon.y(:,:,r_point,plot_c(i)));
            subplot(3,3,i)
            h=pcolor(y_LatLon.lon*90/pi,y_LatLon.lat*90/pi,to_plot);
            set(h, 'EdgeColor', 'none');
            colorbar
            colormap(cmocean('balance',20))
            c_lim=max(abs(to_plot(:)));
            caxis([-c_lim c_lim])
            title(tit{i},'interpreter','latex')
            set(gca,'Box','on');
            set(gca,'fontsize', 18); 
        end
    case 'test_strain_stress' %%%------------ Test Stress Strain 
        
        fig=figure;
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.9, 0.9]);        
        p = uipanel('Parent',fig,'BorderType','none'); 
        if isempty(plot_title)==0
        p.Title = plot_title; 
        p.TitlePosition = 'centertop'; 
        p.FontSize = 18;
        p.FontWeight = 'bold';
        end
        
        div=y_LatLon.y(:,:,r_point,14)+y_LatLon.y(:,:,r_point,18)+y_LatLon.y(:,:,r_point,22);
        tit={'$100\Re(\sigma_{rr}-\lambda\epsilon_{kk}-2\mu\epsilon_{rr})/\textrm{max}(\sigma_{rr})$',...
            '$100\Re(\sigma_{r\theta}-2\mu\epsilon_{r\theta})/\textrm{max}(\sigma_{r\theta})$',...
            '$100\Re(\sigma_{r\phi}-2\mu\epsilon_{r\phi})/\textrm{max}(\sigma_{r\phi})$',...
            '$100\Re(\sigma_{\theta r}-2\mu\epsilon_{\theta r})/\textrm{max}(\sigma_{\theta r})$',...
            '$100\Re(\sigma_{\theta\theta}-\lambda\epsilon_{kk}-2\mu\epsilon_{\theta\theta})/\textrm{max}(\sigma_{\theta\theta})$',...
            '$100\Re(\sigma_{\theta \phi}-2\mu\epsilon_{\theta \phi})/\textrm{max}(\sigma_{\theta \phi})$',...
            '$100\Re(\sigma_{\phi r}-2\mu\epsilon_{\phi r})/\textrm{max}(\sigma_{\phi r})$',...
            '$100\Re(\sigma_{\phi \theta}-2\mu\epsilon_{\phi \theta})/\textrm{max}(\sigma_{\phi \theta})$',...
            '$100\Re(\sigma_{\phi \phi }-\lambda\epsilon_{kk}-2\mu\epsilon_{\phi \phi })/\textrm{max}(\sigma_{\phi \phi })$'};
%          tit={'$100\Re(\sigma_{rr}-\lambda\epsilon_{kk}-2\mu\epsilon_{rr})/\sigma_{rr}$',...
%             '$100\Re(\sigma_{r\theta}-2\mu\epsilon_{r\theta})/\sigma_{r\theta}$',...
%             '$100\Re(\sigma_{r\phi}-2\mu\epsilon_{r\phi})/\sigma_{r\phi}$',...
%             '$100\Re(\sigma_{\theta r}-2\mu\epsilon_{\theta r})/\sigma_{\theta r}$',...
%             '$100\Re(\sigma_{\theta\theta}-\lambda\epsilon_{kk}-2\mu\epsilon_{\theta\theta})/\sigma_{\theta\theta}$',...
%             '$100\Re(\sigma_{\theta \phi}-2\mu\epsilon_{\theta \phi})/\sigma_{\theta \phi}$',...
%             '$100\Re(\sigma_{\phi r}-2\mu\epsilon_{\phi r})/\sigma_{\phi r}$',...
%             '$100\Re(\sigma_{\phi \theta}-2\mu\epsilon_{\phi \theta})/\sigma_{\phi \theta}$',...
%             '$100\Re(\sigma_{\phi \phi }-\lambda\epsilon_{kk}-2\mu\epsilon_{\phi \phi })/\sigma_{\phi \phi }$'};
        for i=1:9
            subplot(3,3,i,'Parent',p)
            if i==1 || i==5 || i==9
                %to_plot=real((2*y_LatLon.mu.*y_LatLon.y(:,:,r_point,13+i)+(Interior_Model.Ks/Interior_Model.mu00R-2/3*y_LatLon.mu).*div))-real(y_LatLon.y(:,:,r_point,4+i));
                to_plot=real((2*y_LatLon.mu.*y_LatLon.y(:,:,r_point,13+i)+(Interior_Model.Ks-2/3*y_LatLon.mu).*div))-real(y_LatLon.y(:,:,r_point,4+i));
                aux=real(y_LatLon.y(:,:,r_point,4+i));
                if max(abs(aux(:)))>0
                    to_plot=to_plot/max(abs(aux(:)));
%                     to_plot=to_plot./aux;
                end
            else
                to_plot=(real(2*y_LatLon.mu.*y_LatLon.y(:,:,r_point,13+i))-real(y_LatLon.y(:,:,r_point,4+i)));
                aux=real(y_LatLon.y(:,:,r_point,4+i));
                if max(abs(aux(:)))>0
                    to_plot=to_plot/max(abs(aux(:)));
%                     to_plot=to_plot./aux;
                end
            end
            h=pcolor(y_LatLon.lon*90/pi,y_LatLon.lat*90/pi,100*to_plot);
            set(h, 'EdgeColor', 'none');
            colorbar
            colormap(cmocean('balance',20))
            c_lim=100*max(abs(to_plot(:)));
            if abs(c_lim) >= 1e-10
                caxis([-c_lim c_lim])
            end
            title(tit{i},'interpreter','latex')
            set(gca,'fontsize', 17);
        end
    otherwise
        error('Wrong Field')
end
% save plot 
set(gcf,'color','w');
if (save_plot)==1
        %export_fig(fig,save_name,'-pdf','-opengl')
        save_name=['Plots_Out/' save_name];
        export_fig(fig,save_name ,'-png','-opengl','-r300')
end 

end

