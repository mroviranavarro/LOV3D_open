%Author: Marc Rovira-Navarro
%title: LatLon_SPH
%Creation: 03/01/2017
%Last Modified: 01/02/2017
%Institution: DTU Space

%Description:
%function used to convert from lat lon data to surface spherical harmonics
%method used in Blais et al Optimization of Spherical Harmonic Trasform
%Computations
%Variables:
    %Input:
        %zonlat 
    %Output: 
        %stokes: stoles coefficinets
            %stokes.clm
            %stokes.slm
            %stokes.lmax: max degree

function [stokes] = LatLon_SPH(zlonlat,varargin)
stokes.lmax=zlonlat.lmax;
lon2=zlonlat.lon;
lat2=zlonlat.lat;
z2=zlonlat.z;
%% COMPUTE POINTS NEEDED FOR inversion
%points in which it needs to be evaluated
N=zlonlat.lmax;
stokes.lmax=N;
%interpolate data if needed goes here. Grid should be as in SPH_LatLon
%% Obtain spherical harmonics
%the harmonic coefficients are obtained by using Driscoll and Healy 1994
%compute legendre polynomials
colat2=90-lat2;
colat2=flip(colat2);
t=cosd(colat2(:,1)')';
z_aux=flipud(z2);
P_lm= Legendre(N,t);
clm=zeros(N+1,N+1);
slm=zeros(N+1,N+1);
qj=zeros(1,length(colat2(:,1)));

for j=0:N
    a_j=0;
    for h=0:(N-1)
        a_j=a_j+1/(2*h+1)*sin((2*h+1)*(j+0.5)*pi/(2*N));
    end
    a_j=1/(N)*sin(pi*(j+0.5)/(2*N))*a_j;
  qj(j+1)=a_j;
end

qj(N+1:2*N)=flip(qj(1:N));
aa2=zeros(2*N+1,N+1);
bb2=zeros(2*N+1,N+1);
%obtain um and vm
for j=0:2*N-1
    for m=0:N
        y=fft([z_aux(j+1,2*N+1:4*N) z_aux(j+1,1:2*N)]);
        aa2(j+1,:)=real(y(1:N+1));
        bb2(j+1,:)=-imag(y(1:N+1));
    end
end

for l=0:N %up to degree N
    for m=0:l
        for j=0:2*N-1 %loop colat
            a_j=qj(j+1);
            LP=P_lm(l+1,m+1,j+1);
            clm(l+1,m+1)=clm(l+1,m+1)+a_j*LP*aa2(j+1,m+1);
            slm(l+1,m+1)=slm(l+1,m+1)+a_j*LP*bb2(j+1,m+1);
        end
        % (-1)^m added bc I am using P_l^m and not P_{l,m}
        clm(l+1,m+1)=(-1)^m*(1/N)/4*clm(l+1,m+1);
        slm(l+1,m+1)=(-1)^m*(1/N)/4*slm(l+1,m+1);
    end
end

stokes.clm=clm; 
stokes.slm=slm; 



end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEEDED FUNCTIONS


