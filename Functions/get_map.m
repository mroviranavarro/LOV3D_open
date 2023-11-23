%% GET_STRESS_STRAIN_MAP
% AUTHOR: M. Rovira-Navarro 
% USE: obtain potential,displacement,stress,strain in lat lon map 
%% INPUT 
    % y_spectra: radial functions
        % y_Spectra.nf: degree of the forcing
        % y_Spectra.mf: order of the forcing 
        % y_Spectra.n: degree of the solution
        % y_Spectra.m: order of the solution 
        % y.y(radial_point,X,mode) 
            % y(radial_point,1,mode): r radial position
            % y(radial_point,2,mode): U radial displacement
            % y(radial_point,3,mode): V tangential displacement
            % y(radial_point,4,mode): W toroidal displacement
            % y(radial_point,5,mode): R normal stress
            % y(radial_point,6,mode): S tangential stress
            % y(radial_point,7,mode): T toroidal stress 
            % y(radial_point,8,mode): \phi gravitational potential
            % y(radial_point,9,mode): \dot\phi potential stress       
            % y(radial_point,10,mode): u_{n,n-1}
            % y(radial_point,11,mode): u_{n,n}
            % y(radial_point,12,mode): u_{n,n+1}
            % y(radial_point,13,mode): \sigma_{n,n,0}       
            % y(radial_point,14,mode): \sigma_{n,n-2,2}
            % y(radial_point,15,mode): \sigma_{n,n-1,2}
            % y(radial_point,16,mode): \sigma_{n,n,2}
            % y(radial_point,17,mode): \sigma_{n,n+1,2}
            % y(radial_point,18,mode): \sigma_{n,n+2,2}      
            % y(radial_point,19,mode): \epsilon_{n,n,0}
            % y(radial_point,20,mode): \epsilon_{n,n-2,2}
            % y(radial_point,21,mode): \epsilon_{n,n-1,2}
            % y(radial_point,22,mode): \epsilon_{n,n,2}
            % y(radial_point,23,mode): \epsilon{n,n+1,2}
            % y(radial_point,24,mode): \epsilon_{n,n+2,2}
%% OUTPUT 
    % y_LatLon: solution in map
        % y_LatLon.nf: degree of the forcing
        % y_LatLon.mf: order of the forcing 
        % y_LatLon.lon: longitude
        % y_LatLon.lat: latitude 
        % y_LatLon.r: radial point
        % y_LatLon.mu(longitude,latitude): shear modulus 
        % y_LatLon.forcing(longitude,latitude): forcing
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
   
%% FUNCTION ------------------------------------------------------------------ 
function [y_LatLon] = get_map(y_Spectra,Interior_Model,varargin)
%% GET SPHERICAL HARMONIC & DERIVATIVES
[SPH]=get_SPH_functions(500,max(y_Spectra.n));
%% INITIALIZE
y_LatLon.nf=y_Spectra.nf; 
y_LatLon.mf=y_Spectra.mf; 
y_LatLon.lon=SPH(1).lon; 
y_LatLon.lat=SPH(1).lat;
y_LatLon.r=y_Spectra.y(:,1,1);
y_LatLon.y=zeros(size(SPH(1).lon,1),size(SPH(1).lon,2),length(y_LatLon.r),22);
%% CONVERT TO MAP 
n_SPH=[];
m_SPH=[];
for i=1:length(SPH)
    n_SPH=[n_SPH SPH(i).n];
    m_SPH=[m_SPH SPH(i).m];
end
%get mu map
mu_variable=Interior_Model.rheology_variable(:,[1 2 4]);
mu_map=zeros(size(SPH(1).lon,1),size(SPH(1).lon,2));
if Interior_Model.mu_variable(3)~=0 || Interior_Model.K_variable(3)~=0  || Interior_Model.eta_variable(3)~=0
    for i=1:size(mu_variable,1)
        ind=find(mu_variable(i,1)==n_SPH & mu_variable(i,2)==m_SPH);
        mu_map=mu_map+mu_variable(i,3)*SPH(ind).Y;
    end
end

% figure; 
% subplot(2,2,1); pcolor(SPH(1).lon,SPH(1).lat,real(SPH(4).Y)) ; colorbar; shading interp
% subplot(2,2,2); pcolor(SPH(1).lon,SPH(1).lat,imag(SPH(4).Y)); colorbar; shading interp
% subplot(2,2,3); pcolor(SPH(1).lon,SPH(1).lat,real(SPH(7).Y)); colorbar; shading interp
% subplot(2,2,4); pcolor(SPH(1).lon,SPH(1).lat,imag(SPH(7).Y)); colorbar; shading interp
% 
% th = linspace(0,pi,50);    % inclination
% phi = linspace(-pi,pi,50); % azimuth
% [th,phi] = meshgrid(th,phi);
% figure; 
% subplot(2,2,1); pcolor(phi,pi/2-th,real(harmonicY(1,1,th,phi))); colorbar; shading interp
% subplot(2,2,2); pcolor(phi,pi/2-th,imag(harmonicY(1,1,th,phi))); colorbar; shading interp
% subplot(2,2,3); pcolor(phi,pi/2-th,real(harmonicY(2,0,th,phi))); colorbar; shading interp
% subplot(2,2,4); pcolor(phi,pi/2-th,imag(harmonicY(2,0,th,phi))); colorbar; shading interp
y_LatLon.mu=(Interior_Model.muC+mu_map);
%y_LatLon.mu=1+mu_map; 
% get forcing map 
y_LatLon.forcing=zeros(size(SPH(1).lon,1),size(SPH(1).lon,2));
ind=find(y_LatLon.nf==n_SPH & y_LatLon.mf==m_SPH);
y_LatLon.forcing=SPH(ind).Y;

for i=1:length(y_Spectra.n)
    n=y_Spectra.n(i);
    m=y_Spectra.m(i);
    ind=find(n==n_SPH & m==m_SPH);
    Y=SPH(ind).Y;
    E=SPH(ind).E;
    G=SPH(ind).G;
    F=SPH(ind).F;
    H=SPH(ind).H;
    [TR1, TR2]=transformation_matrices(n);
    Phi_S=y_Spectra.y(:,8,i);
    U_S=transpose(y_Spectra.y(:,10:12,i));
    Stress_S=transpose(y_Spectra.y(:,13:18,i));
    Strain_S=transpose(y_Spectra.y(:,19:end,i));
    % potential to map 
    y_LatLon.y(:,:,:,1)=y_LatLon.y(:,:,:,1)+reshape(Phi_S,[1,1,length(Phi_S)]).*Y;
    % displacement to map 
    Y_Displacement=TR1(:,:,1)*U_S;
    E_Displacement=TR1(:,:,2)*U_S;
    G_Displacement=TR1(:,:,3)*U_S;
    % Stress & Strain to Map 
    Y_Stress=TR2(:,:,1)*Stress_S;
    E_Stress=TR2(:,:,2)*Stress_S;
    F_Stress=TR2(:,:,3)*Stress_S;
    G_Stress=TR2(:,:,4)*Stress_S;
    H_Stress=TR2(:,:,5)*Stress_S;
    Y_Strain=TR2(:,:,1)*Strain_S;
    E_Strain=TR2(:,:,2)*Strain_S;
    F_Strain=TR2(:,:,3)*Strain_S;
    G_Strain=TR2(:,:,4)*Strain_S;
    H_Strain=TR2(:,:,5)*Strain_S;
    for k=1:size(Stress_S,2)
        for j=1:9
            y_LatLon.y(:,:,k,4+j)=y_LatLon.y(:,:,k,4+j)+...
                            Y_Stress(j,k)*Y+...
                            E_Stress(j,k)*E+...
                            F_Stress(j,k)*F+...
                            G_Stress(j,k)*G+...
                            H_Stress(j,k)*H; 
           y_LatLon.y(:,:,k,13+j)=y_LatLon.y(:,:,k,13+j)+...
                            Y_Strain(j,k)*Y+...
                            E_Strain(j,k)*E+...
                            F_Strain(j,k)*F+...
                            G_Strain(j,k)*G+...
                            H_Strain(j,k)*H;              
        end
        for j=1:3
            y_LatLon.y(:,:,k,1+j)=y_LatLon.y(:,:,k,1+j)+...
                            Y_Displacement(j,k)*Y+...
                            E_Displacement(j,k)*E+...
                            G_Displacement(j,k)*G;
        end
    end
end
end
%% SPHERICAL HARMONICS AND DERIVATIVES
function [SPH]=get_SPH_functions(N,nmax)
Nth=N;
Nlon=N;
%used for numerical differentiation
TH_C=linspace(0,pi,Nth+1);
PHI_C=linspace(-pi,pi,Nlon+1);
ind_i=2:2:Nth;
ind_i2=4:4:Nth-1;
ind_i3=2:2:length(ind_i)-1;
TH_C2=TH_C(ind_i);
PHI_C2=PHI_C(ind_i);
TH=TH_C(ind_i2);
PHI=PHI_C(ind_i2);
Delta_TH1=TH_C(3)-TH_C(1);
Delta_TH2=TH_C(6)-TH_C(2);
Delta_PHI1=PHI_C(3)-PHI_C(1);
Deltha_PHI2=PHI_C(6)-PHI_C(2);
%get Legendre functions
Plm=Legendre(nmax,cos(TH_C)');
i=1;
for n=0:nmax
    for m=-n:n
        n_v(i)=n; 
        m_v(i)=m; 
        i=i+1;
    end
end
NSph=length(n_v);
[phi,th] = meshgrid(PHI,TH);
i=1;
for n=0:nmax
    for m=-n:n
        longi=exp(1i*m*PHI_C);
        SPH(i).lat=pi/2-th;
        SPH(i).lon=phi;
        SPH(i).n=n; 
        SPH(i).m=m; 
        if m==0
            Y=repmat(squeeze(Plm(n+1,abs(m)+1,:)),[1 length(longi)]).*repmat(longi,[length(TH_C) 1]);
        elseif m>0
            Y=repmat((-1)^m*squeeze(Plm(n+1,abs(m)+1,:)),[1 length(longi)]).*repmat(longi,[length(TH_C) 1])/sqrt(2);
        else
            Y=repmat(squeeze(Plm(n+1,abs(m)+1,:)),[1 length(longi)]).*repmat(longi,[length(TH_C) 1])/sqrt(2);
        end          
        
        E=(Y(3:2:end,2:2:end)-Y(1:2:end-1,2:2:end))/Delta_TH1;
        G=(Y(2:2:end,3:2:end)-Y(2:2:end,1:2:end-1))/Delta_PHI1./repmat(sin(TH_C2)',[1 length(PHI_C2)]);
        F=(E(3:2:end,2:2:end-1)-E(1:2:end-2,2:2:end-1))/Delta_TH2;
        H=(G(3:2:end,2:2:end-1)-G(1:2:end-2,2:2:end-1))/Delta_TH2;        
        SPH(i).Y=Y(ind_i2,ind_i2);
        SPH(i).E=E(ind_i3,ind_i3);
        SPH(i).G=G(ind_i3,ind_i3);
        SPH(i).F=F;
        SPH(i).H=H;
        i=i+1;
    end
end
end
%% CHANGE OF BASE MATRICES
function [TR1, TR4]=transformation_matrices(n)
    TR1=zeros(3,3,3);
    TR2=zeros(9,9,5);
    if n>0
    %% RANKR 1 TENSOR
    %----- Y_{n,n-1}
    ind1=1;
    norm=sqrt(n*(2*n+1));
    %e_r
    ind2=1;
    TR1(ind2,ind1,1)=n/norm;% Y
    TR1(ind2,ind1,2)=0;% E
    TR1(ind2,ind1,3)=0;% G
    %e_theta
    ind2=2;
    TR1(ind2,ind1,1)=0;% Y
    TR1(ind2,ind1,2)=1/norm;% E
    TR1(ind2,ind1,3)=0;% G
    %e_phi
    ind2=3;
     % Y
    TR1(ind2,ind1,1)=0;% Y
    TR1(ind2,ind1,2)=0;% E
    TR1(ind2,ind1,3)=1/norm;% G        
    %----- Y_{n,n}
    ind1=2;
    norm=sqrt(n*(n+1));
    %e_r
    ind2=1;
    TR1(ind2,ind1,1)=0;% Y
    TR1(ind2,ind1,2)=0;% E
    TR1(ind2,ind1,3)=0;% G
    %e_theta
    ind2=2;
    TR1(ind2,ind1,1)=0;% Y
    TR1(ind2,ind1,2)=0;% E
    TR1(ind2,ind1,3)=1i/norm;% G
    %e_phi
    ind2=3;
     % Y
    TR1(ind2,ind1,1)=0;% Y
    TR1(ind2,ind1,2)=-1i/norm;% E
    TR1(ind2,ind1,3)=0;% G 
    
    %----- Y_{n,n+1}
    ind1=3;
    norm=sqrt((n+1)*(2*n+1));
    %e_r
    ind2=1;
    TR1(ind2,ind1,1)=-(n+1)/norm;% Y
    TR1(ind2,ind1,2)=0;% E
    TR1(ind2,ind1,3)=0;% G
    %e_theta
    ind2=2;
    TR1(ind2,ind1,1)=0;% Y
    TR1(ind2,ind1,2)=1/norm;% E
    TR1(ind2,ind1,3)=0;% G
    %e_phi
    ind2=3;
     % Y
    TR1(ind2,ind1,1)=0;% Y
    TR1(ind2,ind1,2)=0;% E
    TR1(ind2,ind1,3)=1/norm;% G 
    else
    %----- Y_{n,n+1}
    ind1=3;
    norm=sqrt((n+1)*(2*n+1));
    %e_r
    ind2=1;
    TR1(ind2,ind1,1)=-(n+1)/norm;% Y
    end
    %% RANK 2 TENSOR
    % Zerilli to Y 
    n1_Y=[n-1 n-1 n-1 n n n n+1 n+1 n+1];
    n2_Y=[n-2 n-1 n n-1 n n+1 n n+1 n+2];
    l_Z=[0 2 2 2 2 2];
    n2_Z=[n n-2 n-1 n n+1 n+2];
    TR3=zeros(9,6);
    for n1=n-1:n+1
        for n2=n1-1:n1+1
            ind=find(n1==n1_Y & n2==n2_Y);
            indl0=find(0==l_Z & n2==n2_Z);
            indl2=find(2==l_Z & n2==n2_Z);
            TR3(ind,indl0)=(-1)^(n+n2)*sqrt((2*n1+1))*Wigner6j(1,0,1,n,n1,n2);
            TR3(ind,indl2)=(-1)^(n+n2)*sqrt(5*(2*n1+1))*Wigner6j(1,2,1,n,n1,n2);
        end
    end
    if n>0
    %------------ Y_{n,n-1,n-2} (1)
    ind1=1; 
    if n>1
        norm=sqrt(n*(n-1)*(2*n-1)*(2*n+1));
    else
        norm=1e13;
    end
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=n*(n-1)/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=(n-1)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=(n-1)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=(n-1)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=n/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=(n-1)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=-n^2/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=-1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    
    %------------ Y_{n,n-1,n-1} (2)
    ind1=2; 
    if n>1
        norm=1i*n*sqrt((n-1)*(2*n+1));
    else
        norm=1e13;
    end
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=(1-n)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=-1/norm; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=n^2/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=(n-1)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=n/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    
    %------------ Y_{n,n-1,n} (3)
    ind1=3; 
    norm=n*sqrt(4*n^2-1);
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=-n^2/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-n/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=-n/norm; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=(n-1)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=n/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=(n-1)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=-n^2/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=-1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    
    %------------ Y_{n,n,n-1} (4)
    ind1=4; 
    norm=1i*n*sqrt((2*n+1)*(n+1));
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=-n/norm; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=n/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=1/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=-1/norm; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-1/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=n*(n+1)/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    
    %------------ Y_{n,n,n} (5)
    ind1=5; 
    norm=n*(n+1);
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-1/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=n*(n+1)/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=-1/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=-1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    
    %------------ Y_{n,n,n+1} (6)
    ind1=6; 
    norm=1i*(n+1)*sqrt(n*(2*n+1));
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=(n+1)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-(n+1)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=1/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=-1/norm; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-1/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=n*(n+1)/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    
    %------------ Y_{n,n+1,n} (7)
    ind1=7; 
    norm=(n+1)*sqrt((2*n+3)*(2*n+1));
    %e_r e_r (2)
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=-(n+1)^2/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=(n+1)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=(n+1)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-(n+2)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=-(n+1)/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=-(n+2)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=-(n+1)^2/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=-1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    
    %------------ Y_{n,n+1,n+1} (8)
    ind1=8; 
    norm=1i*(n+1)*sqrt((n+2)*(2*n+1));
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=(n+2)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=-1/norm; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=(n+1)^2/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-(n+2)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=-(n+1)/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    
    %------------ Y_{n,n+1,n+2} (9)
    ind1=9; 
    norm=sqrt((n+1)*(n+2)*(2*n+1)*(2*n+3));
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=(n+1)*(n+2)/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_theta (2)
    ind2=2;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-(n+2)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    %e_r e_phi (3)
    ind2=3;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=-(n+2)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    %e_theta e_r (4)
    ind2=4;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=-(n+2)/norm; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=-(n+1)/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    % e_phi e_r (7)
    ind2=7;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=-(n+2)/norm; % G
    TR2(ind2,ind1,5)=0; % H
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=0; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=0; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=1/norm; % H
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=-(n+1)^2/norm; % Y
    TR2(ind2,ind1,2)=0; % E
    TR2(ind2,ind1,3)=-1/norm; % F
    TR2(ind2,ind1,4)=0; % G
    TR2(ind2,ind1,5)=0; % H
    else 
        %% n==0
    %------------ Y_{n,n+1,n} (7)
    ind1=7; 
    norm=(n+1)*sqrt((2*n+3)*(2*n+1));
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=-(n+1)^2/norm; % Y
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=-(n+1)/norm; % Y
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=-(n+1)^2/norm; % Y
    %------------ Y_{n,n+1,n+1} (8)
    ind1=8; 
    norm=1i*(n+1)*sqrt((n+2)*(2*n+1));
    % e_theta e_phi (6)
    ind2=6;
    TR2(ind2,ind1,1)=(n+1)^2/norm; % Y
    % e_phi e_theta
    ind2=8;
    TR2(ind2,ind1,1)=-(n+1)/norm; % Y
    %------------ Y_{n,n+1,n+2} (9)
    ind1=9; 
    norm=sqrt((n+1)*(n+2)*(2*n+1)*(2*n+3));
    %e_r e_r (1)
    ind2=1;
    TR2(ind2,ind1,1)=(n+1)*(n+2)/norm; % Y
    % e_theta e_theta (5)
    ind2=5;
    TR2(ind2,ind1,1)=-(n+1)/norm; % Y
    %e_phi e_phi (9)
    ind2=9;
    TR2(ind2,ind1,1)=-(n+1)^2/norm; % Y
    end
    TR4(:,:,1)=TR2(:,:,1)*TR3; 
    TR4(:,:,2)=TR2(:,:,2)*TR3; 
    TR4(:,:,3)=TR2(:,:,3)*TR3; 
    TR4(:,:,4)=TR2(:,:,4)*TR3; 
    TR4(:,:,5)=TR2(:,:,5)*TR3; 
end
% test with some plots
% test_ind=9;
% fig=figure;
% set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.9]);
% subplot(3,2,1)
% h=pcolor(SPH(test_ind).lon,SPH(test_ind).lat,real(SPH(test_ind).Y));
% set(h, 'EdgeColor', 'none');
% colorbar
% colormap(cmocean('diff',20))
% c_lim=max(abs(real(SPH(test_ind).Y(:))));
% caxis([-c_lim c_lim])
% title(['Y_{' num2str(SPH(test_ind).n) num2str(SPH(test_ind).m) '}' ],'interpreter','tex')
% subplot(3,2,3)
% h=pcolor(SPH(test_ind).lon,SPH(test_ind).lat,real(SPH(test_ind).E));
% set(h, 'EdgeColor', 'none');
% colorbar
% colormap(cmocean('diff',20))
% c_lim=max(abs(real(SPH(test_ind).E(:))));
% caxis([-c_lim c_lim])
% title(['E_{' num2str(SPH(test_ind).n) num2str(SPH(test_ind).m) '}' ],'interpreter','tex')
% subplot(3,2,4)
% h=pcolor(SPH(test_ind).lon,SPH(test_ind).lat,real(SPH(test_ind).F));
% set(h, 'EdgeColor', 'none');
% colorbar
% colormap(cmocean('diff',20))
% c_lim=max(abs(real(SPH(test_ind).F(:))));
% caxis([-c_lim c_lim])
% title(['F_{' num2str(SPH(test_ind).n) num2str(SPH(test_ind).m) '}' ],'interpreter','tex')
% subplot(3,2,5)
% h=pcolor(SPH(test_ind).lon,SPH(test_ind).lat,real(SPH(test_ind).G));
% set(h, 'EdgeColor', 'none');
% colorbar
% colormap(cmocean('diff',20))
% c_lim=max(abs(real(SPH(test_ind).G(:))));
% caxis([-c_lim c_lim])
% title(['G_{' num2str(SPH(test_ind).n) num2str(SPH(test_ind).m) '}' ],'interpreter','tex')
% subplot(3,2,6)
% h=pcolor(SPH(test_ind).lon,SPH(test_ind).lat,real(SPH(test_ind).H));
% set(h, 'EdgeColor', 'none');
% colorbar
% colormap(cmocean('diff',20))
% c_lim=max(abs(real(SPH(test_ind).H(:))));
% caxis([-c_lim c_lim])
% title(['H_{' num2str(SPH(test_ind).n) num2str(SPH(test_ind).m) '}' ],'interpreter','tex')

