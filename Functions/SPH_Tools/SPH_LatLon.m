%Author: Marc Rovira-Navarro
%title: SPH_LatLon
%Creation: 07/02/2017
%Last Modified: 07/02/2017
%Institution: DTU Space

%Description:
%function used to convert spherical harmonics to lon, lat grid 
%Computations
%Variables:
    %Input:
        %stokes: only in stokes.lmax stokes.clm and stokes.slm needed, it
        %can also be a vector of stokes structures
            %stokes.lmax: maximum harmonic degree
            %stokes.clm: clm coefficnets stokes.clm(l+1,m+1)
            %stokes.slm: slm coefficients stokes.slm(l+1,m+1)
           
            
         %varargin{1}:
                %limits [lon_min lon_max lat_min lat_max]
                    %if not given all map is computed
                %delta_lat 
                    %if not given pi/2*l_max
                %delta_lon 
                    %if not given pi/2*l_max               
    %Output: 
        %zlonlat
            %zlonlat.lon: [lonmin delta_lon lonmax]
            %zlonlat.lat: [latmin delta_lat latmax];
            %zlonlat.z:  pxq grid where p stands for latitude and q for
            %zlonlat.lmax
        
%Needed functions: 
    %- Legendre must be in the same directory 
function [zlonlat] = SPH_LatLon(stokes,varargin)
zlonlat.lmax=stokes.lmax;
%% variables
%default variables 
delta_lon=180/(2*stokes.lmax); %nominal as delta_lon
delta_lat=180/(2*stokes.lmax); 
lon0=-180+delta_lon/2; 
lonM=180-delta_lon/2;
lat0=-90+delta_lat/2; 
latM=90-delta_lon/2;
for i=1:2:length(varargin)
    switch varargin{i}
        case 'limits'
            lon0=varargin{i+1}(1);
            lonM=varargin{i+1}(2);
            lat0=varargin{i+1}(3);
            latM=varargin{i+1}(4);
        case 'delta_lat'
            delta_lat=varargin{i+1};
    end
end

%% define grid
lon=lon0:delta_lon:lonM;
lat=lat0:delta_lat:latM; 
colat=90-lat;

%% compute 
l_max=stokes.lmax;
colat=colat*pi/180;
lon=lon*pi/180;
P_lm=Legendre(l_max,cos(colat)');
gv=zeros(length(colat),length(lon));
for l=1:l_max+1 %loop in dedgree
     for m=1:l %loop in order
         % (-1) added for consistency with definition in my paper
         % % (-1)^m added bc I am using P_l^m and not P_{l,m}
         v1=(-1)^(m-1)*stokes.clm(l,m)*squeeze(P_lm(l,m,:));
         v2=(-1)^(m-1)*stokes.slm(l,m)*squeeze(P_lm(l,m,:));
         for i=1:length(v1)
             gv(i,:)=gv(i,:)+(v1(i)*cos((m-1)*lon)+v2(i)*sin((m-1)*lon));
         end
     end
end
lon=lon*180/pi; 
lon0=min(lon);
lonM=max(lon);
lat0=min(lat);
latM=max(lat);
[lon_g,lat_g] = meshgrid(lon0:delta_lon:lonM,lat0:delta_lat:latM);
zlonlat.lon=lon_g;
zlonlat.lat=lat_g;
zlonlat.z=gv; 

end