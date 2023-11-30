%% SCRIPT USED TO TEST GET_LOVE
%close all
clear all
clc
set(0,'defaulttextInterpreter','latex') 
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);
cd('..')
addpath(genpath(pwd))
Startall = tic;
%% CONSTANTS
G=6.67e-11;
%% ORBITAL PARAMETRES
% We consider Io in its eccentric orbit, note that the tidal amplitude is
% normalized by (omegaR)^2e
omega0=4.1086E-05; %orbital frequency  
T=2*pi/omega0;
Forcing(1).Td=T;
Forcing(1).n=2; 
Forcing(1).m=0; 
Forcing(1).F=3/4*sqrt(1/5); 
Forcing(2).Td=2*pi/omega0;
Forcing(2).n=2; 
Forcing(2).m=-2; 
Forcing(2).F=-7/8*sqrt(6/5);
Forcing(3).Td=2*pi/omega0;
Forcing(3).n=2; 
Forcing(3).m=2; 
Forcing(3).F=1/8*sqrt(6/5);
%% INTERIOR MODEL 
% The interior model is defined in the structure Interior_Model. 
% The interior model can either be specified using dimensional or non-dimensional parameters. 
% In case dimensional parameters are employed, they are non-dimensionalised in get_Love L161
% In this example we use the non-dimensional parameters of a one-layer Io (see Section 2.2 and Table 1)
mu_eff=5.2; %effective shear modulus 
r_ratio=0.53; %ratio between surface and core radius R/R_2
rho_r=1.59; %ratio between core and mantle density 
%Interior_Model.Delta_rho can be use to specify the density contrast between the uppermost
%layer and the layer immediately below. For icy moons Delta_rho=0 can be specified(ocean
%and water have the same density). If it is not specified it simply follows
%from rho_r
Ks_nd=200/60; %ratio between bulk and shear modulus
Interior_Model.MaxTime=3.42; %non-dimensional Maxwell-time 
Interior_Model.eta=Interior_Model.MaxTime/(2*pi); %non-dimensional viscosity, here from the Maxwell time
% uncomment if the model is elastic
% Interior_Model.eta=NaN;
% Interior_Model.MaxTime=NaN;
% here the interior model is build 
%nondimensional model 
Interior_Model.eta0=Interior_Model.eta; 
rho_av=rho_r*r_ratio^3+(1-r_ratio^3);
Interior_Model.R=[r_ratio 1];
Interior_Model.rho=[rho_r 1];
Interior_Model.Ks=Ks_nd; 
Interior_Model.mu=1;
Interior_Model.Gg=3/(4*pi)/mu_eff/rho_av^2; %non-dimensional gravitational constatant
%dimensional parameters
Interior_Model.R0=[r_ratio 1];
Interior_Model.rho0=[rho_r 1];
Interior_Model.Ks0=Ks_nd; 
Interior_Model.mu0=1;
% LATERAL VARIATIONS
% Lateral variations are given in terms of real spherical harmonics.
% In this case, we provide the amplitude of each harmonics in % of peak-to-peak variations with respect the mean properities
nR=2; %degree of lateral variations
mR=0; %order of lateral variations 
variable_mu_p=0; % shear modulus variations
variable_eta_p=50; % viscosity variations 
variable_K_p=0; % bulk modulus variations 
%% NUMERICS
Numerics.Nr=500; %number of radial points 
Numerics.perturbation_order=2; %maximum order to which couplings are considered
Numerics.rheology_cutoff=2; % terms of the rheology that are considerd in the expansion, only relevant for viscoelastic, terms with log10(mu_n^m)-log10(mu_n^m(leading))>=-Numerics.rheology_cutoff are included 
Numerics.load_couplings=2; % rheology couplings are loaded from file that contains all loading coefficients. 
% The file is located in Files_Couplings. If the file does not exist, it is obtained and stored.
Numerics.Nenergy=12; %maximum degree to which the energy dissipation is expanded 
Numerics.parallel=0; %should this be run in parallel?
%% convert lateral variations to SPH
% this part converts the variations given in real spherical harmonics and
% percentage into complex spherical harmonics.
% it can be commented out if the amplitude of the complex spherical harmonics is already given in: Interior_Model.mu_variable, Interior_Model.eta_variable
l_max=30;
Ynm_stokes.clm=zeros(2*l_max,2*l_max);
Ynm_stokes.slm=zeros(2*l_max,2*l_max);
Ynm_stokes.lmax=2*l_max-1;
Ynm_stokes.clm(nR+1,mR+1)=1;
[Ynm_zlonlat] = SPH_LatLon(Ynm_stokes);
Ynm_z=Ynm_zlonlat.z;
Delta=max(Ynm_z(:))-min(Ynm_z(:));
if mR==0
    % shear modulus 
    Interior_Model.mu_variable(1,1)=nR;
    Interior_Model.mu_variable(1,2)=mR;
    Interior_Model.mu_variable(1,3)=variable_mu_p/100/Delta;
    % bulk modulus
    Interior_Model.K_variable(1,1)=nR;
    Interior_Model.K_variable(1,2)=mR;
    Interior_Model.K_variable(1,3)=variable_K_p/100/Delta;
    %viscosity 
    Interior_Model.eta_variable(1,1)=nR;
    Interior_Model.eta_variable(1,2)=mR;
    Interior_Model.eta_variable(1,3)=variable_eta_p/100/Delta;
elseif mR>0
    % shear modulus 
    Interior_Model.mu_variable(1,1)=nR;
    Interior_Model.mu_variable(2,1)=nR;
    Interior_Model.mu_variable(1,2)=mR;
    Interior_Model.mu_variable(2,2)=-mR;    
    Interior_Model.mu_variable(1,3)=sqrt(2)/2*variable_mu_p/100/Delta;
    Interior_Model.mu_variable(2,3)=(-1)^mR*sqrt(2)/2*variable_mu_p/100/Delta;
    %bulk modulus
    Interior_Model.K_variable(1,1)=nR;
    Interior_Model.K_variable(2,1)=nR;
    Interior_Model.K_variable(1,2)=mR;
    Interior_Model.K_variable(2,2)=-mR;    
    Interior_Model.K_variable(1,3)=sqrt(2)/2*variable_K_p/100/Delta;
    Interior_Model.K_variable(2,3)=(-1)^mR*sqrt(2)/2*variable_K_p/100/Delta;
    % viscosity 
    Interior_Model.eta_variable(1,1)=nR;
    Interior_Model.eta_variable(2,1)=nR;
    Interior_Model.eta_variable(1,2)=mR;
    Interior_Model.eta_variable(2,2)=-mR;    
    Interior_Model.eta_variable(1,3)=sqrt(2)/2*variable_eta_p/100/Delta;
    Interior_Model.eta_variable(2,3)=(-1)^mR*sqrt(2)/2*variable_eta_p/100/Delta;
else
    % shear modulus 
    Interior_Model.mu_variable(1,1)=nR;
    Interior_Model.mu_variable(2,1)=nR;
    Interior_Model.mu_variable(1,2)=-mR;
    Interior_Model.mu_variable(2,2)=mR;    
    Interior_Model.mu_variable(1,3)=-1i*(-1)^mR*sqrt(2)/2*variable_mu_p/100/Delta;
    Interior_Model.mu_variable(2,3)=1i*sqrt(2)/2*variable_mu_p/100/Delta;
    % bulk modulus
    Interior_Model.K_variable(1,1)=nR;
    Interior_Model.K_variable(2,1)=nR;
    Interior_Model.K_variable(1,2)=-mR;
    Interior_Model.K_variable(2,2)=mR;    
    Interior_Model.K_variable(1,3)=-1i*(-1)^mR*sqrt(2)/2*variable_K_p/100/Delta;
    Interior_Model.K_variable(2,3)=1i*sqrt(2)/2*variable_K_p/100/Delta;
    %viscosity
    Interior_Model.eta_variable(1,1)=nR;
    Interior_Model.eta_variable(2,1)=nR;
    Interior_Model.eta_variable(1,2)=-mR;
    Interior_Model.eta_variable(2,2)=mR;    
    Interior_Model.eta_variable(1,3)=-1i*(-1)^mR*sqrt(2)/2*variable_eta_p/100/Delta;
    Interior_Model.eta_variable(2,3)=1i*sqrt(2)/2*variable_eta_p/100/Delta;
end
% prepares a spherically-symmetric model for comparison   
Interior_Model_Uni=Interior_Model;
Interior_Model_Uni.eta_variable(:,3)=0; 
Interior_Model_Uni.mu_variable(:,3)=0; 
Interior_Model_Uni.K_variable(:,3)=0; 
%% GET LOVE NUMBERS
%obtains the Love number spectra
if Numerics.parallel==1
    parfor i=1:length(Forcing)
        % spherically-symmetric model 
        [Love_Spectra_Uni(i),y_Uni(i),Interior_Model_UniU(i)]=get_Love(Interior_Model_Uni,Forcing(i),Numerics,'verbose','out_file','T1');
        % lateral variations
        if i==1 %plot rheology maps
            [Love_Spectra(i),y(i),Interior_ModelU(i)]=get_Love(Interior_Model,Forcing(i),Numerics,'plot_rheology_spectra',1,'verbose','out_file','T1');
        else
            [Love_Spectra(i),y(i),Interior_ModelU(i)]=get_Love(Interior_Model,Forcing(i),Numerics,'plot_rheology_spectra',0,'verbose','out_file','T1');
        end
    end
else
    for i=1:length(Forcing)
        % spherically-symmetric model 
        [Love_Spectra_Uni(i),y_Uni(i),Interior_Model_UniU(i)]=get_Love(Interior_Model_Uni,Forcing(i),Numerics,'verbose','out_file','T1');
        % lateral variations
        if i==1 %plot rheology maps
            [Love_Spectra(i),y(i),Interior_ModelU(i)]=get_Love(Interior_Model,Forcing(i),Numerics,'plot_rheology_spectra',1,'verbose','out_file','T1');
        else
            [Love_Spectra(i),y(i),Interior_ModelU(i)]=get_Love(Interior_Model,Forcing(i),Numerics,'plot_rheology_spectra',0,'verbose','out_file','T1');
        end
    end
end
%% GET ENERGY SPECTRA
% obtains the energy dissipation spectra
% can be commentted if elastic
[Energy_Spectra]=get_energy(y,Numerics,Forcing,Interior_ModelU(end),'verbose',1,'out_file','T1');
[Energy_Spectra_Uni]=get_energy(y_Uni,Numerics,Forcing,Interior_Model_UniU(end),'verbose',1,'out_file','T1');
Energy_Spectra_Norm=Energy_Spectra_Uni; 
Energy_Spectra_Norm.energy_integral_v=Energy_Spectra_Norm.energy_integral_v/Energy_Spectra_Uni.energy_integral_v(1); 
Energy_Spectra_Diff=Energy_Spectra;
Energy_Spectra_Diff.energy_integral_v=(Energy_Spectra.energy_integral_v-Energy_Spectra_Uni.energy_integral_v)/Energy_Spectra_Uni.energy_integral_v(1);
disp(['Time spent on all calculations ' num2str(toc(Startall)) ' s'])
%% PLOT ENERGY SPECTRA, 
% can be commented if elastic 
% total tidal heating for unifrom model 
plot_energy_map(Energy_Spectra_Norm,'type','total','projection','mollweide','title','Uniform Normalized');
% total tidal heating for the model with lateral variations
plot_energy_map(Energy_Spectra_Norm,'type','total','projection','mollweide','title','Variable');
% difference in tidal heating between the models with and without lateral
% variations
plot_energy_map(Energy_Spectra_Diff,'type','difference','projection','cut','title','Difference','limits',[-0.4 0.4],'label','$\Delta\dot{e}/\dot{e}_0^u$');
%% GET MAPS AND PLOT
% Here the y functions given in spherical harmonics are mapped to a lat lon and some relavnt tensors are plotted.  
% (2,0) Component ------------------
[y_LatLon] = get_map(y(1),Interior_ModelU(end));
[y_LatLonUni] = get_map(y_Uni(1),Interior_Model_UniU(end));
yDif=y_LatLon; 
yDif.y=y_LatLon.y-y_LatLonUni.y;
%plot the gravitational potential, displacements and strain and stress tensors for the uniform model 
plot_map(y_LatLonUni,Interior_Model_UniU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Uniform Model (2,0)') 
%plot the gravitational potential, displacements and strain and stress tensors for the variable model 
plot_map(y_LatLon,Interior_ModelU(end),'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Variable Model (2,0)');
%plot the difference in gravitational potential, displacements and strain and stress tensors between the variable and the uniform models. 
plot_map(yDif,Interior_ModelU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Difference (2,0)')
% these two plots can be used to check that the stress and strain fields are correctly computed, the residuals should be close to 0. 
% If not, a higher perturbation order or Numerics.rheology_cutoff should be used 
plot_map(y_LatLon,Interior_ModelU(end),'field_name','test_strain_stress','radial_point',floor(Numerics.Nr*0.5),'plot_title','(2,0)');
plot_map(y_LatLonUni,Interior_Model_UniU(end),'field_name','test_strain_stress','radial_point',floor(Numerics.Nr*0.5),'plot_title','(2,0) Uniform');
% (2,2) Component -------------------------
[y_LatLon] = get_map(y(2),Interior_ModelU(end));
[y_LatLonUni] = get_map(y_Uni(2),Interior_Model_UniU(end));
yDif=y_LatLon; 
yDif.y=y_LatLon.y-y_LatLonUni.y;
%plot the gravitational potential, displacements and strain and stress tensors for the uniform model 
plot_map(y_LatLonUni,Interior_Model_UniU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Uniform Model (2,2)') 
%plot the gravitational potential, displacements and strain and stress tensors for the variable model 
plot_map(y_LatLon,Interior_ModelU(end),'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Variable Model (2,2)');
%plot the difference in gravitational potential, displacements and strain and stress tensors between the variable and the uniform models. 
plot_map(yDif,Interior_ModelU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Difference (2,2)')
% these two plots can be used to check that the stress and strain fields are correctly computed, the residuals should be close to 0. 
% If not, a higher perturbation order or Numerics.rheology_cutoff should be used 
plot_map(y_LatLon,Interior_ModelU(end),'field_name','test_strain_stress','radial_point',floor(Numerics.Nr*0.5),'plot_title','(2,2)');
plot_map(y_LatLonUni,Interior_Model_UniU(end),'field_name','test_strain_stress','radial_point',floor(Numerics.Nr*0.5),'plot_title','(2,2) Uniform');
%% TEST ENERGY DISSIPATION 
% compute energy dissipation using Eq. (34) and  (38)
E_k_Uni=0;
E_k=0;
for i=1:length(Forcing)
    n_f=Forcing(i).n;
    m_f=Forcing(i).m;
    for j=1:length(Forcing)
        ind1=find(Love_Spectra_Uni(j).n==n_f & Love_Spectra_Uni(j).m==m_f);
        ind2=find(Love_Spectra(j).n==n_f & Love_Spectra(j).m==m_f);
        if isempty(ind1)==0
            E_k_Uni=E_k_Uni-Forcing(i).F*Forcing(j).F*imag(Love_Spectra_Uni(j).k(ind1)); 
        end
        if isempty(ind2)==0
            E_k=E_k-Forcing(i).F*Forcing(j).F*imag(Love_Spectra(j).k(ind2)); 
        end
    end
end
% 
e_01_Uni=2*pi*10/(Interior_Model_UniU(1).Gg)*E_k_Uni;
e_02_Uni=4*pi*Energy_Spectra_Uni.energy_integral(1);
disp([ 'Difference between energy computation with method 1 and 2 for the uniform model ' num2str((e_01_Uni-e_02_Uni)/e_01_Uni*100) ' \%'])
%variable
e_01=2*pi*10/(Interior_ModelU(1).Gg)*E_k;
e_02=4*pi*Energy_Spectra.energy_integral(1);
disp([ 'Difference between energy computation with method 1 and 2 for a model with lateral variations ' num2str((e_01-e_02)/e_01*100) ' \%'])
%% TRANSFORM DIMENSIONAL UNITS
% an example showing how to transform to dimensional units for Io
G=6.67e-11;
ecc=4e-3; %Io eccentricity 
rho_0=3244; %density of uppermost layer
R_0=1821.3e3; %surface radius 
mu_0=60e9; %shear modulus 
scale=omega0^2*R_0^2*ecc*rho_0/mu_0; %scale for the tidal potential 
%uniform 
E_total_Uni=R_0*scale^2*mu_0/T*Energy_Spectra_Uni.energy_integral(1); %(in W/m2)
%variable 
E_total=R_0*scale^2*mu_0/T*Energy_Spectra.energy_integral(1); %(in W/m2)
disp([ 'Tidal heating spherically-symmetric model ' num2str(E_total_Uni) ' W/m^2'])
disp([ 'Tidal heating laterally-heterogenous model ' num2str(E_total) ' W/m^2'])