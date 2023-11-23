%% PLOT_RHEOLOGY_MAP
% function used to plot map of rheology pattern
%% INPUT

%% OUTPUT 

%%
function [varargout] = plot_rheology_map(plot_rheology_spectra,mu_zlonlat,eta_zlonlat,MaxTime_zlonlat, ...
                                         Interior_Model,Cmu_zlonlat,m_g,n_g, ...
                                         muI_aux,muI_SPH,muR_aux,muR_SPH, ...
                                         mu00I,non_zero_indexes,varargin)
if plot_rheology_spectra == 1
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.9]);
    %shear modulus
    ax(1)=subplot(4,2,1)
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(mu_zlonlat.lon,mu_zlonlat.lat,1+mu_zlonlat.z); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    colormap(ax(1),cmocean('gray'))
    set(0,'defaulttextInterpreter','none')
    m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',12,'gridcolor',[0.8 0.8 0.8],'linewidth',2,'linestyle','--','gridcolor','k');
    cb=colorbar;
    set(0,'defaulttextInterpreter','latex')
    cb.Label.String='$\mu/\mu_{0,0}$';
    cb.Label.Interpreter = 'latex';
    set(gca,'fontsize', 18);
    set(gcf,'color','w');
    %viscosity
    ax(2)=subplot(4,2,2)
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(mu_zlonlat.lon,mu_zlonlat.lat,1+eta_zlonlat.z); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    colormap(ax(2),cmocean('matter'))
    set(0,'defaulttextInterpreter','none')
    m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',12,'gridcolor',[0.8 0.8 0.8],'linewidth',2,'linestyle','--','gridcolor','k');
    cb=colorbar;
    set(0,'defaulttextInterpreter','latex')
    cb.Label.String='$\eta/\eta_{0,0}$';
    cb.Label.Interpreter = 'latex';
    set(gca,'fontsize', 18);
    set(gcf,'color','w');
    % Maxwell time non-dimensional
    ax(3)=subplot(4,2,3)
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(mu_zlonlat.lon,mu_zlonlat.lat,MaxTime_zlonlat.z); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    colormap(ax(3),cmocean('matter'))
    set(0,'defaulttextInterpreter','none')
    m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',12,'gridcolor',[0.8 0.8 0.8],'linewidth',2,'linestyle','--','gridcolor','k');
    cb=colorbar;
    set(0,'defaulttextInterpreter','latex')
    cb.Label.String='$\frac{\tau_M}{\eta_0 / T \mu_0 }$';
    cb.Label.Interpreter = 'latex';
    set(gca,'fontsize', 18);
    set(gcf,'color','w');
    % Maxwell time dimensional
    ax(3)=subplot(4,2,4)
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(mu_zlonlat.lon,mu_zlonlat.lat,MaxTime_zlonlat.z*Interior_Model.MaxTime); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    colormap(ax(3),cmocean('matter'))
    set(0,'defaulttextInterpreter','none')
    m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',12,'gridcolor',[0.8 0.8 0.8],'linewidth',2,'linestyle','--','gridcolor','k');
    cb=colorbar;
    set(0,'defaulttextInterpreter','latex')
    cb.Label.String='$\frac{\tau_M}{T}$';
    cb.Label.Interpreter = 'latex';
    set(gca,'fontsize', 18);
    set(gcf,'color','w');
    %Real part  of mu
    ax(4)=subplot(4,2,5)
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(mu_zlonlat.lon,mu_zlonlat.lat,real(Cmu_zlonlat.z)); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    colormap(ax(4),cmocean('gray'))
    set(0,'defaulttextInterpreter','none')
    m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',12,'gridcolor',[0.8 0.8 0.8],'linewidth',2,'linestyle','--','gridcolor','k');
    cb=colorbar;
    set(0,'defaulttextInterpreter','latex')
    cb.Label.String='$Re(\hat\mu) / \mu_{0}$';
    cb.Label.Interpreter = 'latex';
    set(gca,'fontsize', 18);
    set(gcf,'color','w');
    ax(5)=subplot(4,2,6)
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(mu_zlonlat.lon,mu_zlonlat.lat,imag(Cmu_zlonlat.z)); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    colormap(ax(5),flip(cmocean('matter')))
    set(0,'defaulttextInterpreter','none')
    m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',12,'gridcolor',[0.8 0.8 0.8],'linewidth',2,'linestyle','--','gridcolor','k');
    cb=colorbar;
    set(0,'defaulttextInterpreter','latex')
    cb.Label.String='$Im(\hat\mu)/ \mu_{0}$';
    cb.Label.Interpreter = 'latex';
    set(gca,'fontsize', 18);
    set(gcf,'color','w');
    % spectra real part 
    ax(6)=subplot(4,2,7)
    s=pcolor(m_g-0.5,n_g-0.5,log10(abs(muR_aux/muR_SPH.clm(1,1))))
    s.EdgeColor = [1 1 1];
    cb=colorbar
    c1=(cmocean('gray',13));
    c2=(cmocean('matter',6));
    colorsM=[c1(1:7,:); flip(c2(1:3,:))];
    colormap(ax(6),colorsM)
    cb.Label.String='$\textrm{log}\left(\textrm{Re}(\hat\mu)/\textrm{Re}(\mu)_{0,0}\right)$';
    cb.Label.Interpreter = 'latex';
    caxis([-10 0])
    ylabel('Degree')
    xlabel('Order')
    grid on
    set(gca,'fontsize', 18);
    % spectra imaginary part
    ax(7)=subplot(4,2,8)
    s=pcolor(m_g-0.5,n_g-0.5,log10(abs(muI_aux/muI_SPH.clm(1,1))))
    s.EdgeColor = [1 1 1];
    cb=colorbar
    colormap(ax(7),colorsM)
    caxis([-10 0])
    ylabel('Degree')
    xlabel('Order')
    cb.Label.String='$\textrm{log}\left(\textrm{Im}(\hat\mu)/\textrm{Im}(\mu)_{0,0}\right)$';
    cb.Label.Interpreter = 'latex';
    set(gca,'fontsize', 18);
    
    % % test that mu I is the same
    % th = linspace(0,pi,500);    % inclination
    % phi = linspace(0,2*pi,500); % azimuth
    % [th,phi] = meshgrid(th,phi); 
    % muI_map=zeros(size(th)); 
    % 
    % muI_map=mu00I*harmonicY(0,0,th,phi);
    % for i=1:length(non_zero_indexes)
    % if n_g(non_zero_indexes(i))>0
    %     muI_map=muI_map+muI_aux(non_zero_indexes(i))*harmonicY(n_g(non_zero_indexes(i)),m_g(non_zero_indexes(i)),th,phi);
    % end
    % end
    % aux2=real(muI_map);    
    % fig=figure;
    % set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.9]);
    % m_proj('mollweide','clongitude',0);
    % set(0,'defaulttextInterpreter','none') 
    % obj=m_pcolor(phi*180/pi-180,90-th*180/pi,real(muI_map)); 
    % set(0,'defaulttextInterpreter','latex') 
    % shading interp;
    % colormap(flip(cmocean('matter')))
    % set(0,'defaulttextInterpreter','none')
    % m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',12,'gridcolor',[0.8 0.8 0.8],'linewidth',2,'linestyle','--','gridcolor','k');
    % cb=colorbar;
    % set(0,'defaulttextInterpreter','latex')
    % cb.Label.Interpreter = 'latex';
    % set(gca,'fontsize', 18);
    % set(gcf,'color','w');   

elseif plot_rheology_spectra==2 %make plots to save
    % Map of the Maxwell-Time 
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.4, 0.4]);
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(mu_zlonlat.lon,mu_zlonlat.lat,MaxTime_zlonlat.z*Interior_Model.MaxTime); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    colormap(cmocean('matter'))
    set(0,'defaulttextInterpreter','none')
    m_grid('xaxislocation','middle','xtick',8,'ytick',7,'fontsize',12,'gridcolor',[0.8 0.8 0.8],'linewidth',2,'linestyle','--','gridcolor','k');
    cb=colorbar;
    set(0,'defaulttextInterpreter','latex')
    cb.Label.String='$2\pi\tau_M/T$';
    cb.Label.Interpreter = 'latex';
    set(gca,'fontsize', 24);
    set(gcf,'color','w');
    %save_name=['Test_visco/RheologyMap_R' num2str(variable_eta(1)) num2str(variable_eta(2))];
    %export_fig(fig,save_name,'-pdf','-opengl')    
    % Real part of the complex modulus 
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.6]);
    s=pcolor(m_g-0.5,n_g-0.5,log10(abs(muR_aux/muR_SPH.clm(1,1))))
    s.EdgeColor = [1 1 1];
    c1=(cmocean('gray',104));
    c2=(cmocean('matter',48));
    colorsM=[c1(1:56,:); flip(c2(1:24,:))];
    colormap(colorsM)
    ylabel('Degree')
    xlabel('Order')
    xlim([-10.5 10.5])
    ylim([-0.5 10.5])
    set(gca,'fontsize', 24);
    label_colorbar='$\textrm{Re}(\hat\mu)_n^m/\textrm{Re}(\mu)_{0}^{0}$';
    color_bar_log([-10 0],-10:2:0,label_colorbar);
    %save_name=['Test_visco/RheologyReal_R' num2str(variable_eta(1)) num2str(variable_eta(2))];
    set(gcf,'color','w');
    %export_fig(fig,save_name,'-pdf','-opengl')    
    % Imaginary part of the complex modulus
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.6]);
    s=pcolor(m_g-0.5,n_g-0.5,log10(abs(muI_aux/muI_SPH.clm(1,1))))
    s.EdgeColor = [1 1 1];
    colormap(colorsM)
    ylabel('Degree')
    xlabel('Order')
    xlim([-10.5 10.5])
    ylim([-0.5 10.5])
    set(gca,'fontsize', 24);
    label_colorbar='$\textrm{Im}(\hat\mu)_n^m/\textrm{Im}(\mu)_{0}^{0}$';
    color_bar_log([-10 0],-10:2:0,label_colorbar);
    set(gcf,'color','w');
    %save_name=['Test_visco/RheologyImag_R' num2str(variable_eta(1)) num2str(variable_eta(2))];
    %export_fig(fig,save_name,'-pdf','-opengl')
end

%%check that conversion has been done correctly---------SHOULD BE AT THE
%%END OF THE FIRST IF STATEMENT
% %              muI_SPH2=muI_SPH; 
% %              %muI_SPH2.slm(:)=0;
% %              muR_SPH2=muR_SPH; 
% %              %muR_SPH2.slm(:)=0;
% %     %        muI_SPH2.clm(:)=0;
% %     %        muI_SPH2.clm(1:2,1)=muI_SPH.clm(1:2,1);
% %             [muI_zlonlat2]=SPH_LatLon(muI_SPH2);
% %             [muR_zlonlat2]=SPH_LatLon(muR_SPH2);  
% %            % complex-spherical harmonics to map 
% %             [SPH]=get_SPH_functions(500,l_max);
% %             n_SPH=[];
% %             m_SPH=[];
% %             for i=1:length(SPH)
% %                 n_SPH=[n_SPH SPH(i).n];
% %                 m_SPH=[m_SPH SPH(i).m];
% %             end
% %             mu_variable=Interior_Model.rheology_variable(:,[1 2 4]);
% %             mu_map=zeros(size(SPH(1).lon,1),size(SPH(1).lon,2));
% %             if Interior_Model.mu_variable(3)~=0 || Interior_Model.K_variable(3)~=0  || Interior_Model.eta_variable(3)~=0
% %                 for i=1:size(mu_variable,1)
% %                     ind=find(mu_variable(i,1)==n_SPH & mu_variable(i,2)==m_SPH);
% %                     mu_map=mu_map+mu00R*mu_variable(i,3)*SPH(ind).Y;
% %                 end
% %             end
% % 
% %             lUI=max(muI_zlonlat.z(:));
% %             lDI=min(muI_zlonlat.z(:));
% %             lUR=max(muR_zlonlat.z(:));
% %             lDR=min(muR_zlonlat.z(:)); 
% % 
% %             % imaginary component 
% %             fig=figure;
% %             set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.9]);
% %             %original ---
% %             subplot(4,2,1)
% %             pcolor( muI_zlonlat.lon, muI_zlonlat.lat, muI_zlonlat.z)
% %             title('Original Im')
% %             shading interp;
% %             colorbar
% %             caxis([lDI lUI]);
% %             subplot(4,2,2)
% %             pcolor( muR_zlonlat.lon, muR_zlonlat.lat, muR_zlonlat.z)
% %             title('Original Re')
% %             shading interp;
% %             colorbar
% %             caxis([lDR lUR]);
% %             %converted to sph harm
% %             subplot(4,2,3)
% %             pcolor( muI_zlonlat2.lon, muI_zlonlat2.lat, muI_zlonlat2.z)
% %             shading interp;
% %             title('Converted to spherical harmonics Im')
% %             colorbar
% %             caxis([lDI lUI]);
% %             subplot(4,2,4)
% %             pcolor( muR_zlonlat2.lon, muR_zlonlat2.lat, muR_zlonlat2.z)
% %             shading interp;
% %             title('Converted to spherical harmonics Re')
% %             colorbar
% %             caxis([lDR lUR]);
% %             %Transformed to complex spherical harmonics and cut Im
% %             subplot(4,2,5)
% %             pcolor(SPH(1).lon, SPH(1).lat,imag(mu00R*mu00+mu_map))
% %             title('Transformed to complex spherical harmonics and cut Im')
% %             shading interp;
% %             colorbar
% %             caxis([lDI lUI]);
% %             subplot(4,2,6)
% %             pcolor(SPH(1).lon, SPH(1).lat,real(mu00R*mu00+mu_map))
% %             title('Transformed to complex spherical harmonics and cut Re')
% %             shading interp;
% %             colorbar
% %             caxis([lDR lUR]);
% %             %Difference
% %             subplot(4,2,7)
% %             pcolor( muI_zlonlat2.lon, muI_zlonlat2.lat,(muI_zlonlat2.z-muI_zlonlat.z)./muI_zlonlat.z)
% %             shading interp;
% %             colorbar
% %             subplot(4,2,8)
% %             pcolor( muR_zlonlat2.lon, muR_zlonlat2.lat,(muR_zlonlat2.z-muR_zlonlat.z)./muR_zlonlat.z)
% %             shading interp;
% %             colorbar