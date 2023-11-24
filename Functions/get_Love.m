%% GET_LOVE
% AUTHOR: M. Rovira-Navarro & A. Veenstra
% USE: obtain the tidal love number spectra for an interior model with lateral variations
%% INPUT 
% Interior_Model: Structure containing the interior model information
% Dimensional (indicated with a 0) or nondimensional parameters can be
% provided. If dimensional parameters are provided, they are
% nondiemensionalised inside the code
        %Interior_Model.R0: radius
            %R0(1) core radius (solid+liquid core)
            %R0(2) surface radius
        %Interior_Model.rho0: layers density        
            %rho0(1) density of interior layer (solid+liquid core)
            %rho0(2) density of outermost solid layer
        %Interior_Model.Delta_rho0: Density difference between the liquid core and overlying solid layer. If not provided it is computed assuming that the two innermost layers have rho0(1).
        %Interior_Model.mu0: shear modulus of the the outermost layer
        %Interior_Model.Ks0: bulk modulus of the outermost layer
        %Interior_Model.eta0: viscosity of the outermost layer
        % Interio_Model.mu_variable: shear modulus variations
            %mu_variable(:,1): degree of variation 
            %mu_variable(:,2): order of variation 
            %mu_variable(:,3): amplitude of the variation (mu_l^m/mu^0_0)
        % Interio_Model.K_variable: bulk modulus variations 
            %K_variable(:,1): degree of variation 
            %K_variable(:,2): order of variation 
            %K_variable(:,3):  amplitude of bulk modulus variations (K_l^m/K^0_0)
        % Interio_Model.eta_variable: viscosity
            %eta_variable(:,1): degree of variation 
            %eta_variable(:,2): order of variation 
            %eta_variable(:,3):  amplitude of viscosity variations (eta_l^m/eta^0_0)    
        % Interior_Model.rheology_variable:  complex shear and bulkd modulus lateral variations (assigned inside the code)
            %rheology_variable(:,1): degree of variation 
            %rheology_variable(:,2): order of variation 
            %rheology_variable(:,3):  amplitude of bulk modulus variations (K_l^m/K^0_0), not used or tested yet 
            %rheology_variable(:,4):  amplitude of complex shear modulus variations (mu_l^m/mu^0_0) 
        % Interior_Model.muC: complex shear modulus (assigned inside the code)   

    % Forcing: Structure containing forcing information
        % Forcing.Td: forcing period
        % Forcing.n: degree of the forcing 
        % Forcing.m: order of the forcing 
        % Forcing.F: amplitude of the component 

    % Numerics: Structure containing the numerical information
        % Numerics.Nr: number of radial points
        % Numerics.perturbation_order: maximum perturbation order. Default 2
        % Numerics.rheology_cutoff: determines which terms of the rheology are included (only relevant for viscoelastic).
            % terms with log10(mu_n^m)-log10(mu_n^m(leading))>=-Numerics.rheology_cutoff are included 
            % Default 0 (only leading terms)
        % Numerics.load_couplings: 
            % (0) compute coupling coefficients from scratch
            % (1) load coupling coefficients from file that contains  exactly the same coupling coefficintes. If the file does not exist, it is generated and stored in Files_Coupling
            % (2) load coupling coefficintes from a file that contains all couplung coefficients a given degree. If the file does not exist, it is generated and stored in Files_Coupling 
        % Numerics.Nrheo_max: Maximum degree lateral variations for which coupling coefficients are computed. Only used if Numerics.load_couplings=2.

    % optional variables 
        %plot_rheology_spectra: plot rheology spectra (1)
        %verbose: print information in screen (1)
        %out_file: print output in file named out_file
%% OUTPUT 
    % Love_Spectra: Love number spectra
        % Love_Spectra.nf: degree of the forcing 
        % Love_Spectra.mf: order of the forcing 
        % Love_Spectra.n: degree of the solution
        % Love_Spectra.m: order of the solution
        % Love_Spectra.k: gravity Love numbers
        % Love_Spectra.h: radial displacement Love numbers
    % y_rad: radial functions
        % y_rad.nf: degree of the forcing
        % y_rad.mf: order of the forcing 
        % y_rad.n: degree of the solution
        % y_rad.m: order of the solution 
        % y_rad.y(radial_point,X,mode) 
            % y_rad.y(radial_point,1,mode): r radial position
            % y_rad.y(radial_point,2,mode): U radial displacement
            % y_rad.y(radial_point,3,mode): V tangential displacement
            % y_rad.y(radial_point,4,mode): R normal stress
            % y_rad.y(radial_point,5,mode): S tangential stress
            % y_rad.y(radial_point,6,mode): \phi gravitational potential
            % y_rad.y(radial_point,7,mode): \dot\phi potential stress
            % y_rad.y(radial_point,8,mode): W toroidal displacement
            % y_rad.y(radial_point,9,mode): T toroidal stress        
            % y_rad.y(radial_point,10,mode): u_{n,n-1}
            % y_rad.y(radial_point,11,mode): u_{n,n}
            % y_rad.y(radial_point,12,mode): u_{n,n+1}
            % y_rad.y(radial_point,13,mode): \sigma_{n,n,0}       
            % y_rad.y(radial_point,14,mode): \sigma_{n,n-2,2}
            % y_rad.y(radial_point,15,mode): \sigma_{n,n-1,2}
            % y_rad.y(radial_point,16,mode): \sigma_{n,n,2}
            % y_rad.y(radial_point,17,mode): \sigma_{n,n+1,2}
            % y_rad.y(radial_point,18,mode): \sigma_{n,n+2,2}      
            % y_rad.y(radial_point,19,mode): \epsilon_{n,n,0}
            % y_rad.y(radial_point,20,mode): \epsilon_{n,n-2,2}
            % y_rad.y(radial_point,21,mode): \epsilon_{n,n-1,2}
            % y_rad.y(radial_point,22,mode): \epsilon_{n,n,2}
            % y_rad.y(radial_point,23,mode): \epsilon{n,n+1,2}
            % y_rad.y(radial_point,24,mode): \epsilon_{n,n+2,2}
   % varargout{1}: Interior_Model updated with nondimensionalized parameters       
%% FUNCTION ------------------------------------------------------------------   
function [Love_Spectra,y_rad,varargout] = get_Love(Interior_Model,Forcing,Numerics,varargin)
Interior_Model.Gg0=6.67428E-11;
Interior_Model.Gg0=6.67e-11;
uniform=1; 
elastic=0;
plot_rheology_spectra=0; 
verbose=0; 
out_file=0; 
%% (0) Optional Inputs Defaults -----------------------------------------------------
if isfield(Numerics,'Nr')==0
    Numerics.Nr=100;
end
if isfield(Numerics,'rheology_cutoff')==0
    Numerics.rheology_cutoff=0;
end
if isfield(Numerics,'load_couplings')==0
    Numerics.load_couplings=1;
end
if isfield(Numerics,'perturbation_order')==0
    Numerics.perturbation_orders=2;
end
if isfield(Interior_Model,'eta0')==0 && isfield(Interior_Model,'eta')==0
        elastic=1; 
        Interior_Model.eta0=NaN;
end
if isfield(Interior_Model,'eta_variable')==1 
    if max(abs(Interior_Model.eta_variable(:,3)))>0
        uniform=0;
    end
end
if isfield(Interior_Model,'mu_variable')==1
    if max(abs(Interior_Model.mu_variable(:,3)))>0
        uniform=0;
    end
end
if isfield(Interior_Model,'K_variable')==1
    if max(abs(Interior_Model.eta_variable(:,3)))>0
        uniform=0;
    end
end
for k = 1:length(varargin)
    if strcmpi(varargin{k},'plot_rheology_spectra')
        plot_rheology_spectra=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'verbose')
        verbose=1; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'out_file')
        out_file=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end

%% (0) Nondimensional ------------------------------------------------------------------
if verbose==1
    disp('##############################')
    disp('Strating run...');
    time_start=tic;
end

if isfield(Interior_Model,'Gg')==0
    Interior_Model.R=Interior_Model.R0/Interior_Model.R0(end);
    Interior_Model.rho=Interior_Model.rho0/Interior_Model.rho0(end);
    Interior_Model.Ks=Interior_Model.Ks0/Interior_Model.mu0(end);
    Interior_Model.mu=Interior_Model.mu0/Interior_Model.mu0(end);
    Interior_Model.Gg=Interior_Model.Gg0*Interior_Model.rho0(end)^2*Interior_Model.R0(end)^2/Interior_Model.mu0(end);
    Interior_Model.eta=Interior_Model.eta0/(Interior_Model.mu0*Forcing.Td);
    Interior_Model.MaxTime=2*pi*Interior_Model.eta0/(Forcing.Td*Interior_Model.mu0);
end

if isnan(Interior_Model.eta0)==1 || isnan(Interior_Model.eta)==1
    elastic=1; 
end

% Compute average density, required to get h and l Love numbers
rho_av=4/3*pi*Interior_Model.R0(1)^3*Interior_Model.rho0(1); 
for i=2:length(Interior_Model.R)
    rho_av=rho_av+4/3*pi*Interior_Model.rho0(i)*(Interior_Model.R0(i)^3-Interior_Model.R0(i-1)^3);
end
Interior_Model.rho_av0=rho_av/(4/3*pi*Interior_Model.R0(end)^3);
Interior_Model.rho_av=Interior_Model.rho_av0/Interior_Model.rho0(end);
Interior_Model.gs0=4/3*pi*Interior_Model.Gg0*Interior_Model.rho_av0*Interior_Model.R0(end);
Interior_Model.gs=4/3*pi*Interior_Model.Gg*Interior_Model.rho_av*Interior_Model.R(end);

if isfield(Interior_Model,'Delta_rho')==0
    Interior_Model.Delta_rho=Interior_Model.rho(1)-Interior_Model.rho(2);
end
%% (1) GET RHEOLOGY ------------------------------------------------------------------
l_max=25; %the higher the better the transfomration but also slower 
if elastic==0 % viscoelastic model
    if uniform==0 %lateral variations 
        % get shear modulus in lateral, longitude grid 
        mu_stokes.clm=zeros(2*l_max,2*l_max);
        mu_stokes.slm=zeros(2*l_max,2*l_max);
        mu_stokes.lmax=2*l_max-1;
        for i=1:size(Interior_Model.mu_variable,1)
            % Convert to real spherical harmonics
            if Interior_Model.mu_variable(i,2)>0
                nR=Interior_Model.mu_variable(i,1);
                mR=Interior_Model.mu_variable(i,2);
                amplitude=2/sqrt(2)*Interior_Model.mu_variable(i,3);
                mu_stokes.clm(nR+1,mR+1)=amplitude;
            elseif Interior_Model.mu_variable(i,2)<0               
            else
                nR=Interior_Model.mu_variable(i,1);
                mR=Interior_Model.mu_variable(i,2);
                amplitude=Interior_Model.mu_variable(i,3);
                mu_stokes.clm(nR+1,mR+1)=amplitude;
            end     
        end
        [mu_zlonlat] = SPH_LatLon(mu_stokes);
        % get viscosity in lateral, longitude grid
        eta_stokes.clm=zeros(2*l_max,2*l_max);
        eta_stokes.slm=zeros(2*l_max,2*l_max);
        eta_stokes.lmax=2*l_max-1;
        for i=1:size(Interior_Model.eta_variable,1)
            % Convert to real spherical harmonics
            if Interior_Model.eta_variable(i,2)>0
                nR=Interior_Model.eta_variable(i,1);
                mR=Interior_Model.eta_variable(i,2);
                amplitude=2/sqrt(2)*Interior_Model.eta_variable(i,3);
                eta_stokes.clm(nR+1,mR+1)=amplitude;
            elseif Interior_Model.eta_variable(i,2)<0               
            else
                nR=Interior_Model.eta_variable(i,1);
                mR=Interior_Model.eta_variable(i,2);
                amplitude=Interior_Model.eta_variable(i,3);
                eta_stokes.clm(nR+1,mR+1)=amplitude;
            end
        end
        [eta_zlonlat] = SPH_LatLon(eta_stokes);
        % Compute the Maxwell time 
        MaxTime_zlonlat.lon=mu_zlonlat.lon;
        MaxTime_zlonlat.lat=mu_zlonlat.lat;
        MaxTime_zlonlat.lmax=mu_zlonlat.lmax;
        MaxTime_zlonlat.z=(1+eta_zlonlat.z)./(1+mu_zlonlat.z);
        % Compute the complex Shear modulus
        Cmu_zlonlat.lon=mu_zlonlat.lon;
        Cmu_zlonlat.lat=mu_zlonlat.lat;
        Cmu_zlonlat.lmax=mu_zlonlat.lmax;
        Cmu_zlonlat.z=(1+mu_zlonlat.z)./(1-1i./MaxTime_zlonlat.z/Interior_Model.MaxTime);
        % obtain spherical harmonics coefficients
        % real part 
        muR_zlonlat.lon=mu_zlonlat.lon;
        muR_zlonlat.lat=mu_zlonlat.lat;
        muR_zlonlat.lmax=mu_zlonlat.lmax;
        muR_zlonlat.z=real(Cmu_zlonlat.z);
        [muR_SPH]=LatLon_SPH(muR_zlonlat);
        %imaginary part
        muI_zlonlat.lon=mu_zlonlat.lon;
        muI_zlonlat.lat=mu_zlonlat.lat;
        muI_zlonlat.lmax=mu_zlonlat.lmax;
        muI_zlonlat.z=imag(Cmu_zlonlat.z);
        [muI_SPH]=LatLon_SPH(muI_zlonlat);               
        % set components in the right format for the code (complex spherical harmonics) 
        vec_n=0:1:2*l_max-1;
        vec_m=-2*l_max+1:1:2*l_max-1;
        [m_g, n_g]=meshgrid(vec_m,vec_n);
        muR_aux=zeros(size(m_g));
        muI_aux=zeros(size(m_g));
        for i=1:length(vec_n) %degree 
            for j=1:vec_n(i)+1 %order
                if vec_n(j)==0
                    indexN=find(vec_n(j)==m_g & vec_n(i)==n_g);
                    muR_aux(indexN)=muR_aux(indexN)+muR_SPH.clm(i,j);
                    muI_aux(indexN)=muI_aux(indexN)+muI_SPH.clm(i,j);
                else
                    % Note: Here the slm component is set automatically to 0
                    % This is done because only lateral variations with pattern corresponding to  real spherical harmonics with m>0 have been considered. 
                    % The *0 should be deletated if this changes. 
                    indexN1=find(vec_n(j)==m_g & vec_n(i)==n_g); %positive m
                    indexN2=find(-vec_n(j)==m_g & vec_n(i)==n_g); %negative m
                    % REAL mu
                    %positive m
                    muR_aux(indexN1)=muR_aux(indexN1)+1/sqrt(2)*muR_SPH.clm(i,j);
                    muR_aux(indexN1)=muR_aux(indexN1)-(-1)^vec_n(j)*1i/sqrt(2)*muR_SPH.slm(i,j)*0;
                    %negative m
                    muR_aux(indexN2)=muR_aux(indexN2)+(-1)^vec_n(j)/sqrt(2)*muR_SPH.clm(i,j);
                    muR_aux(indexN2)=muR_aux(indexN2)+1i/sqrt(2)*muR_SPH.slm(i,j)*0;
                    % IMAGINARY mu 
                    muI_aux(indexN1)=muI_aux(indexN1)+1/sqrt(2)*muI_SPH.clm(i,j);
                    muI_aux(indexN1)=muI_aux(indexN1)-(-1)^vec_n(j)*1i/sqrt(2)*muI_SPH.slm(i,j)*0;
                    %negative m
                    muI_aux(indexN2)=muI_aux(indexN2)+(-1)^vec_n(j)/sqrt(2)*muI_SPH.clm(i,j);
                    muI_aux(indexN2)=muI_aux(indexN2)+1i/sqrt(2)*muI_SPH.slm(i,j)*0;
                end
            end
        end
        % 0,0 component
        index0=find(n_g==0 & m_g==0);
        mu00R=muR_aux(index0);
        mu00I=muI_aux(index0);
        mu00=mu00R+1i*mu00I;
        %mu00=1+1i*mu00I/mu00R;
        Interior_Model.muC=mu00;
        %find components that are bigger than some quantity, this is an approximation 
        muI_aux2=log10(abs(muI_aux/muI_SPH.clm(1,1)));
        muR_aux2=log10(abs(muR_aux/muR_SPH.clm(1,1)));
        muI_aux2(index0)=log10(0);
        muR_aux2(index0)=log10(0); 
        muI_max=max(muI_aux2(:));
        muR_max=max(muR_aux2(:));
        mu_max=max([muR_max muI_max]);
        % redefine Gg
        %Interior_Model.Gg=Interior_Model.Gg/mu00R;
        non_zero_indexes=find(muI_aux2-mu_max>=-Numerics.rheology_cutoff | muR_aux2-mu_max>=-Numerics.rheology_cutoff);
        k=1;
        for i=1:length(non_zero_indexes)
            if n_g(non_zero_indexes(i))>0
                Interior_Model.rheology_variable(k,1)=n_g(non_zero_indexes(i)); 
                Interior_Model.rheology_variable(k,2)=m_g(non_zero_indexes(i)); 
                Interior_Model.rheology_variable(k,4)=(muR_aux(non_zero_indexes(i))+1i*muI_aux(non_zero_indexes(i)));
                %Interior_Model.rheology_variable(k,4)=(muR_aux(non_zero_indexes(i))+1i*muI_aux(non_zero_indexes(i)))/mu00R;
                k=k+1;
            end
        end
        Interior_Model.mu00R=mu00R;
        % PLOT THE RHEOLOGY STRUCTURE
        plot_rheology_map(plot_rheology_spectra,mu_zlonlat,eta_zlonlat, ...
                          MaxTime_zlonlat,Interior_Model,Cmu_zlonlat, ...
                          m_g, n_g,muI_aux,muI_SPH,muR_aux,muR_SPH, ...
                          mu00I,non_zero_indexes)

    else
        muC=(1-1i/Interior_Model.MaxTime)^(-1);
        mu00R=real(muC);
        mu00I=imag(muC);
        mu00=mu00R+1i*mu00I;
        %mu00=1+1i*mu00I/mu00R;
        Interior_Model.muC=mu00;
        Interior_Model.mu00R=mu00R;
        % redefine Gg
        %Interior_Model.Gg=Interior_Model.Gg/mu00R;
        Interior_Model.rheology_variable=[0 0 0 0];
    end
else %elastic
    if uniform==1
        Interior_Model.rheology_variable=[0 0 0 0];
        Interior_Model.muC=Interior_Model.mu;
        Interior_Model.mu00R=Interior_Model.mu;
    else
        Interior_Model.rheology_variable(:,[1 2])=Interior_Model.mu_variable(:,[1 2]);
        Interior_Model.rheology_variable(:,3)=0;
        Interior_Model.rheology_variable(:,4)=Interior_Model.mu_variable(:,3);
        Interior_Model.muC=Interior_Model.mu;
        Interior_Model.mu00R=Interior_Model.mu;
    end
end
Interior_Model.lambda=Interior_Model.Ks-2/3*Interior_Model.muC;
%Interior_Model.lambda=Interior_Model.Ks/Interior_Model.mu00R-2/3*Interior_Model.muC;
if uniform==1
    Interior_Model.rheology_variable=[0 0 0 0]; 
end
%% (0.1) File for printing -------------------------------------------------------------------
str_rheo=[];
for i=1:size(Interior_Model.rheology_variable,1)
    if i==size(Interior_Model.rheology_variable,1)
        str_rheo=[str_rheo num2str(Interior_Model.rheology_variable(i,1)) '_' num2str(Interior_Model.rheology_variable(i,2))];
    else
        str_rheo=[str_rheo num2str(Interior_Model.rheology_variable(i,1)) '_' num2str(Interior_Model.rheology_variable(i,2)) '__'];
    end
end
if out_file~=0
    out_file_name=['Files_Out/' out_file '__' str_rheo '_per' num2str(Numerics.rheology_cutoff) '__forc__' num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.txt'];
    FID=fopen(out_file_name,'w');
    out_file_name_reference_y=['Files_Out/y__' str_rheo '_per' num2str(Numerics.rheology_cutoff) '__forc__' num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.mat'];
end
%% (2) OBTAIN COUPLINGS --------------------------------------------------------------
if verbose==1
    disp('Obtaining Couplings....')
    tic
end
if uniform==0
    % Coupling file name    
    if Numerics.load_couplings==2 %files 
        Numerics.Nrheo_max=max(Interior_Model.rheology_variable(:,1));
        coupling_file_name=['Files_Coupling/L__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.mat'];
    else
    coupling_file_name=['Files_Coupling/L__rheo__' str_rheo '__forc__' num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.mat'];
    end

    % Load specific coupling file
    if isfile(coupling_file_name)==1 && abs(Interior_Model.rheology_variable(1,4))>0 && Numerics.load_couplings==1
        if verbose==1
            disp(['Coupling Loaded from: ' coupling_file_name])
        end
        load(coupling_file_name)
    % Check if a coupling file with sufficient modes exist and load it 
    elseif  Numerics.load_couplings==2 
        % look if there is a file that contains the couplings
        coupling_file_name_search=['Files_Coupling/L__Nrheomax__*__forc__' num2str(Forcing.n) '_' num2str(Forcing.m) '_per*.mat'];
        possible_couplings_files=dir(coupling_file_name_search);
        file_found=0;
        i=1; 
        while file_found==0 && i<=length(possible_couplings_files)
            per_start=strfind(possible_couplings_files(i).name,'_per')+4;
            per_end=strfind(possible_couplings_files(i).name,'.mat')-1;
            perturbation_order_file=str2num(possible_couplings_files(i).name(per_start:per_end));
            rheo_start=strfind(possible_couplings_files(i).name,'Nrheomax__')+10;
            rheo_end=strfind(possible_couplings_files(i).name,'__forc')-1;
            Nrheomax_file=str2num(possible_couplings_files(i).name(rheo_start:rheo_end));
            if Nrheomax_file>=Numerics.Nrheo_max && perturbation_order_file>=Numerics.perturbation_order
                file_found=1; 
                coupling_file_name=['Files_Coupling/' possible_couplings_files(i).name];
            end
            i=i+1;
        end
        if verbose==1
            disp(['Coupling Loaded from: ' coupling_file_name])
        end
        if file_found==1 %coupling file exist 
            load(coupling_file_name)
        else %coumpute couplings 
                if verbose==1
                    disp(['File ' coupling_file_name 'not found. Computing all coupling coefficients, this might take some time..'])
                    Couplings=get_couplings_all(Numerics.perturbation_order,Numerics.Nrheo_max,Forcing,'verbose');
                else
                    Couplings=get_couplings_all(Numerics.perturbation_order,Numerics.Nrheo_max,Forcing,'verbose');
                end
            % store coupling file
            save(coupling_file_name,'Couplings')
        end
        % remove terms with 0 amplitude 
        non_zero_rheo=find(abs(Interior_Model.rheology_variable(:,4))>0);
        Interior_Model.rheology_variable=Interior_Model.rheology_variable(non_zero_rheo,:); 
        % retrieve the coupling coefficients that are required
        Couplings=retrieve_couplings(Numerics.perturbation_order,Interior_Model.rheology_variable(:,1:2),Forcing,Couplings);
    else % compute couplings for specific case
        Couplings=get_couplings(Numerics.perturbation_order,Interior_Model.rheology_variable(:,1:2),Forcing);
        save(coupling_file_name,'Couplings')
        if verbose==1
            disp(['Time Spent ' num2str(toc) ' s'])
        end
    end
else
    Couplings.n_s=Forcing.n;
    Couplings.m_s=Forcing.m;
    Couplings.order=0;
    Couplings.Coup=zeros(1,1,27,1); 
end
Nsol=length(Couplings.n_s);
%% PRINT MODEL INFORMATION IN SCREEN & FILE
if verbose==1
    disp('------------- INTERIOR MODEL -----------')
    disp('AVERAGE PROPERTIES (dimensional)')
    disp(['Layer#     R[m]    rho[kg.m^{-3}]    mu[Pa]    K[Pa]    eta[Pa.s]' ])
    disp(['1    ' num2str(Interior_Model.R0(1),'%10.5e') '    ' num2str(Interior_Model.rho0(1),'%10.5e') '    0'  '    -'   ])
    for i=1:length(Interior_Model.R)-1
        disp([num2str(i+1) '    ' num2str(Interior_Model.R0(i+1),'%10.5e') '    ' num2str(Interior_Model.rho0(i+1),'%10.5e') '    ' num2str(Interior_Model.mu0(i),'%10.5e')  '    ' num2str(Interior_Model.Ks0(i),'%10.5e') ' ' num2str(Interior_Model.eta0(i),'%10.5e')   ])
    end
    disp('AVERAGE PROPERTIES (non-dimensional)')
    disp(['Layer#     R[-]    rho[-]    mu[-]    K [-]     eta[-]' ])
    disp(['1    ' num2str(Interior_Model.R(1)) '    ' num2str(Interior_Model.rho(1)) '    0'  '    -'   ])
    for i=1:length(Interior_Model.R)-1
        disp([num2str(i+1) '    ' num2str(Interior_Model.R(i+1)) '    ' num2str(Interior_Model.rho(i+1)) '    ' num2str(Interior_Model.mu(i))  '    ' num2str(Interior_Model.Ks(i)) ' ' num2str(Interior_Model.eta(i))   ])
    end
    disp('RHEOLOGY  VARIATIONS')
    disp('Shear Modulus')
    if isfield(Interior_Model,'mu_variable')==1
        disp('(n,m)           amplitude[mu_n^m/mu_0^0]')
        for i=1:size(Interior_Model.mu_variable,1)
            disp(['(' num2str(Interior_Model.mu_variable(i,1)) ',' num2str(Interior_Model.mu_variable(i,2)) ')   '  num2str(Interior_Model.mu_variable(i,3),'%10.5e') ])
        end
    else
        disp('None')
    end 
    
    disp('Bulk Modulus')
    if isfield(Interior_Model,'K_variable')==1
        disp('(n,m)           amplitude[K_n^m/K_0^0]')
        for i=1:size(Interior_Model.K_variable,1)
            disp(['(' num2str(Interior_Model.K_variable(i,1)) ',' num2str(Interior_Model.K_variable(i,2)) ')   '  num2str(Interior_Model.K_variable(i,3),'%10.5e') ])
        end
    else
        disp('None')
    end
    
    disp('Viscosity')
    if isfield(Interior_Model,'eta_variable')==1
        disp('(n,m)          amplitude[eta_n^m/eta_0^0]')
        for i=1:size(Interior_Model.eta_variable,1)
            disp(['(' num2str(Interior_Model.eta_variable(i,1)) ',' num2str(Interior_Model.eta_variable(i,2)) ')   '  num2str(Interior_Model.eta_variable(i,3),'%10.5e') ])
        end
    else
        disp('None')
    end
    
    disp('Complex Shear Modulus')    
    if isfield(Interior_Model,'rheology_variable')==1
        disp('(n,m)         amplitude[\hat\mu_n^m/mu_0^0]')
        for i=1:size(Interior_Model.rheology_variable,1)
            disp(['(' num2str(Interior_Model.rheology_variable(i,1)) ',' num2str(Interior_Model.rheology_variable(i,2)) ')    '  num2str(Interior_Model.rheology_variable(i,4),'%10.5e') ])
        end
    else
        disp('None')
    end
    disp(' ')
    disp('----------- FORCING ----------')
    disp('TIDAL POTENTIAL')
    disp('Period [s]')
    disp(num2str(Forcing.Td,'%10.5e'))
    disp('(n,m)')
    disp(['(' num2str(Forcing.n) ',' num2str(Forcing.m) ')'])
    disp(' ')
    disp('-----------  RESPONSE SPECTRUM ----------- ')
    for j=0:Numerics.perturbation_order
    index_order=find(Couplings.order==j);
    str_h=[];
    for k=1:length(index_order)
        str_h=[str_h '(' num2str(Couplings.n_s(index_order(k))) ',' num2str(Couplings.m_s(index_order(k))) '),   '];
    end
    disp([num2str(j) 'th Order Modes'])
    disp([num2str(length(index_order)) ' modes'])
    disp(str_h)
    end
    disp('----------- NUMERICAL INFORMATION ----------')
    disp(['Number of Modes ' num2str(Nsol)])
    disp(['Radial Points '  num2str(Numerics.Nr)])
end
% FILE
if out_file~=0
    fprintf(FID, '%s\n','------------- INTERIOR MODEL -----------');
    fprintf(FID, '%s\n', 'AVERAGE PROPERTIES (dimensional)');
    fprintf(FID, '%s\n', ['Layer#     R[m]    rho[kg.m^{-3}]    mu[Pa]    K[Pa]    eta[Pa.s]' ]);
    fprintf(FID, '%s\n', ['1    ' num2str(Interior_Model.R0(1),'%10.5e') '    ' num2str(Interior_Model.rho0(1),'%10.5e') '    0'  '    -'   ]);
    for i=1:length(Interior_Model.R)-1
        fprintf(FID, '%s\n', [num2str(i+1) '    ' num2str(Interior_Model.R0(i+1),'%10.5e') '    ' num2str(Interior_Model.rho0(i+1),'%10.5e') '    ' num2str(Interior_Model.mu0(i),'%10.5e')  '    ' num2str(Interior_Model.Ks0(i),'%10.5e') ' ' num2str(Interior_Model.eta0(i),'%10.5e')   ]);
    end
    fprintf(FID, '%s\n', 'AVERAGE PROPERTIES (non-dimensional)');
    fprintf(FID, '%s\n', ['Layer#     R[-]    rho[-]    mu[-]    K [-]     eta[-]' ]);
    fprintf(FID, '%s\n', ['1    ' num2str(Interior_Model.R(1)) '    ' num2str(Interior_Model.rho(1)) '    0'  '    -'   ]);
    for i=1:length(Interior_Model.R)-1
        fprintf(FID, '%s\n', [num2str(i+1) '    ' num2str(Interior_Model.R(i+1)) '    ' num2str(Interior_Model.rho(i+1)) '    ' num2str(Interior_Model.mu(i))  '    ' num2str(Interior_Model.Ks(i)) ' ' num2str(Interior_Model.eta(i))   ]);
    end
    fprintf(FID, '%s\n', 'RHEOLOGY  VARIATIONS');
    fprintf(FID, '%s\n', ['Rheology cutoff =        ' num2str(Numerics.rheology_cutoff)]);
    fprintf(FID, '%s\n', 'Shear Modulus');
    if isfield(Interior_Model,'mu_variable')==1
        fprintf(FID, '%s\n', '(n,m)           amplitude[mu_n^m/mu_0^0]');
        for i=1:size(Interior_Model.mu_variable,1)
            fprintf(FID, '%s\n', ['(' num2str(Interior_Model.mu_variable(i,1)) ',' num2str(Interior_Model.mu_variable(i,2)) ')   '  num2str(Interior_Model.mu_variable(i,3),'%10.5e') ]);
        end
    else
        fprintf(FID, '%s\n', 'None');
    end 
    
    fprintf(FID, '%s\n', 'Bulk Modulus');
    if isfield(Interior_Model,'K_variable')==1
        fprintf(FID, '%s\n', '(n,m)           amplitude[K_n^m/K_0^0]');
        for i=1:size(Interior_Model.K_variable,1)
            fprintf(FID, '%s\n', ['(' num2str(Interior_Model.K_variable(i,1)) ',' num2str(Interior_Model.K_variable(i,2)) ')   '  num2str(Interior_Model.K_variable(i,3),'%10.5e') ]);
        end
    else
        fprintf(FID, '%s\n', 'None')
    end
    
    fprintf(FID, '%s\n', 'Viscosity');
    if isfield(Interior_Model,'eta_variable')==1
        fprintf(FID, '%s\n', '(n,m)          amplitude[eta_n^m/K_0^0]');
        for i=1:size(Interior_Model.eta_variable,1)
            fprintf(FID, '%s\n', ['(' num2str(Interior_Model.eta_variable(i,1)) ',' num2str(Interior_Model.eta_variable(i,2)) ')   '  num2str(Interior_Model.eta_variable(i,3),'%10.5e') ]);
        end
    else
        fprintf(FID, '%s\n', 'None');
    end
    
    fprintf(FID, '%s\n', 'Complex Shear Modulus');
    
    if isfield(Interior_Model,'rheology_variable')==1
        fprintf(FID, '%s\n', '(n,m)         amplitude[\hat\mu_n^m/K_0^0]');
        for i=1:size(Interior_Model.rheology_variable,1)
            fprintf(FID, '%s\n', ['(' num2str(Interior_Model.rheology_variable(i,1)) ',' num2str(Interior_Model.rheology_variable(i,2)) ')    '  num2str(Interior_Model.rheology_variable(i,4),'%10.5e') ]);
        end
    else
        fprintf(FID, '%s\n', 'None');
    end
    fprintf(FID, '%s\n', ' ');
    fprintf(FID, '%s\n', '----------- FORCING ----------');
    fprintf(FID, '%s\n', 'TIDAL POTENTIAL');
    fprintf(FID, '%s\n', 'Period [s]');
    fprintf(FID, '%s\n', num2str(Forcing.Td,'%10.5e'));
    fprintf(FID, '%s\n', '(n,m)');
    fprintf(FID, '%s\n', ['(' num2str(Forcing.n) ',' num2str(Forcing.m) ')']);
    fprintf(FID, '%s\n', ' ');
    fprintf(FID, '%s\n', '-----------  RESPONSE SPECTRUM ----------- ');
    for j=0:Numerics.perturbation_order
    index_order=find(Couplings.order==j);
    str_h=[];
    for k=1:length(index_order)
        str_h=[str_h '(' num2str(Couplings.n_s(index_order(k))) ',' num2str(Couplings.m_s(index_order(k))) '),   '];
    end
    fprintf(FID, '%s\n',[num2str(j) 'th Order Modes']);
    fprintf(FID, '%s\n',[num2str(length(index_order)) ' modes']);
    fprintf(FID, '%s\n',str_h);
    end
    fprintf(FID, '%s\n', '----------- NUMERICAL INFORMATION ----------');
    fprintf(FID, '%s\n', ['Number of Modes ' num2str(Nsol)]);
    fprintf(FID, '%s\n', ['Radial Points '  num2str(Numerics.Nr)]);
end
%% (3) PROPAGATE AND OBTAIN SOLUTION ------------------------------------------------------------
if verbose==1
    tStart = tic;
end
y_sol=get_solution(Interior_Model,Forcing,Numerics,Couplings,verbose,out_file);
if verbose==1
    disp(['Time Spent: ' num2str(toc(tStart)) 's'])
end
%% (4) PRINT RESPONSE IN SCREEN AND FILE 
if verbose==1
    Nmodes=length(Couplings.n_s);
    disp('---------- RESPONSE ----------------')
    disp('k2 Love numbers')
    disp('(n,m)           k_n^m')
    for i=1:Nmodes
        if Couplings.n_s(i)==Forcing.n && Couplings.m_s(i)==Forcing.m
            k2=y_sol(end,8,i)-1;
        else
            k2=y_sol(end,8,i);
        end
        disp(['(' num2str(Couplings.n_s(i)) ',' num2str(Couplings.m_s(i)) ')           '  num2str(k2,'%10.5e')   ])
    end
end
% FILE
if out_file~=0
    Nmodes=length(Couplings.n_s);
    fprintf(FID, '%s\n','#####################################');
    fprintf(FID, '%s\n','---------- RESPONSE ----------------');
    fprintf(FID, '%s\n','k Love numbers');
    fprintf(FID, '%s\n','(n,m)           k_n^m');
    for i=1:Nmodes
        if Couplings.n_s(i)==Forcing.n && Couplings.m_s(i)==Forcing.m
            k2=y_sol(end,8,i)-1;
        else
            k2=y_sol(end,8,i);
        end
        fprintf(FID, '%s\n',['(' num2str(Couplings.n_s(i)) ',' num2str(Couplings.m_s(i)) ')           '  num2str(k2,'%10.5e')   ]);
    end
    save(out_file_name_reference_y,"y_sol");
end
%% (5) REARRANGE THE SOLUTION  
Nmodes=length(Couplings.n_s);
y_rad.y=y_sol; 
Love_Spectra.nf=Forcing.n;
Love_Spectra.mf=Forcing.m;
y_rad.nf=Forcing.n;
y_rad.mf=Forcing.m;
for i=1:Nmodes
    Love_Spectra.n(i)=Couplings.n_s(i);
    Love_Spectra.m(i)=Couplings.m_s(i);
    Love_Spectra.order(i)=Couplings.order(i);
    Love_Spectra.k(i)=y_sol(end,8,i);
    Love_Spectra.h(i)=-Interior_Model.gs*y_sol(end,2,i);
    Love_Spectra.l(i)=-Interior_Model.gs*y_sol(end,3,i);
    y_rad.n(i)=Couplings.n_s(i);
    y_rad.m(i)=Couplings.m_s(i);
end
varargout{1}=Interior_Model;
if verbose==1
    disp(['Run completed in ' num2str(toc(time_start))])
    disp('##############################')
end
end

%% FUNCTIONS -----------
function [Couplings]=retrieve_couplings(max_order,variations,Forcing,Couplings)
% for the given rheology compute the modes that are excited using the
% selection rules
modes=[];
rheo=unique((variations(:,[1 2])),'rows');
order=0; 
%start with the forcing 
mode0(1)=Forcing.n;
mode0(2)=Forcing.m;
mode0(3)=1;
mode0(4)=0; 
mode0(5)=0;  
modes=[modes mode0];
index_old=1;
order=order+1;
while order<max_order+1
    index_new=[]; 
    for k=1:size(rheo,1)
        modeR(1)=rheo(k,1);
        modeR(2)=rheo(k,2);
        for j=1:length(index_old)
            modes_aux=next_coupling(modes(index_old(j),:),modeR,index_old(j),order);
            index_new=[index_new size(modes,1)+(1:1:size(modes_aux,1))];
            modes=[modes; modes_aux];
        end
    end
    index_old=index_new; 
    order=order+1; 
end

modes2=unique(modes(:,[5, 1:3]),'rows');
modes3=unique(modes(:,[1:2]),'rows');
ind_f=find(modes2(2:end,2)==Forcing.n & modes2(2:end,3)==Forcing.m);
if isempty(ind_f)==0
    modes2(1,1)=modes2(ind_f(1)+1,1);
end
modes4=[];
for j=1:length(modes3)
    ind=find(modes2(:,2)==modes3(j,1) & modes2(:,3)==modes3(j,2));
    modes4(j,:)=[modes3(j,:), min(modes2(ind,1))];
end
n_s=modes4(:,1);
m_s=modes4(:,2);
order=modes4(:,3);

modes_indexes=length(n_s); 
rheo_indexes=size(variations,1);
% find the rheology indexes
for i=1:size(variations,1)
    rheo_indexes(i)=find(variations(i,1)==Couplings.n_r & variations(i,2)==Couplings.m_r);
end
% find the modes indexes. 
for i=1:length(n_s)
    modes_indexes(i)=find(n_s(i)==Couplings.n_s & m_s(i)==Couplings.m_s);
end
% select coefficients
Couplings.n_s=Couplings.n_s(modes_indexes);
Couplings.m_s=Couplings.m_s(modes_indexes);
Couplings.order=order; 
Couplings.n_r=variations(:,1);
Couplings.m_r=variations(:,2);
Couplings.Coup=Couplings.Coup(modes_indexes,modes_indexes,:,rheo_indexes);
end
%% find next set of couplings using selection rules 
function modes=next_coupling(mode0,modeR,parent,order)
n0=mode0(1); 
m0=mode0(2);
ST=mode0(3); 
n1=modeR(1); 
m1=modeR(2);
kmode=1; 
for j=0:1/2*(n0+n1-abs(n0-n1))
    m_new=m0+m1;
    n_new=abs(n0-n1)+2*j; 
    if ST==1
       ST_new=1; 
    else
        ST_new=2; 
    end
    if abs(m_new)<n_new+1
    modes(kmode,1)=n_new;    
    modes(kmode,2)=m_new;
    modes(kmode,3)=ST_new;
    modes(kmode,4)=parent;
    modes(kmode,5)=order;
    kmode=kmode+1;
    end
end
for j=0:1/2*(n0+n1-abs(n0-n1))-1
    m_new=m0+m1;
    n_new=abs(n0-n1)+2*j+1;
    if ST==1
        ST_new=2; 
    else
        ST_new=1; 
    end
    if abs(m_new)<n_new+1 && abs(m_new)+abs(m0)+abs(m1)~=0
    modes(kmode,1)=n_new;    
    modes(kmode,2)=m_new;
    modes(kmode,3)=ST_new;
    modes(kmode,4)=parent;
    modes(kmode,5)=order;
    kmode=kmode+1; 
    end
end
end
