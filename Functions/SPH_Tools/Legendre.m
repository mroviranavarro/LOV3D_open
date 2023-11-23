%Author: Marc Rovira-Navarro
%title: Legendre
%Creation: 02/01/2017
%Last Modified: 02/01/2017
%Institution: DTU Space

%Description:
%function used to compute legendre polynomials value in blocks as in Rapp
%1982
%Variables
    %INPUT
        %l_max: max harmonic degree
    
     %OUTPUT
         %P_lm: associated legendre polynomial value
            %P_lm(l,m,t)=P_(l+1,m+1)(t)


function [P_lm] = Legendre(l_max,t)
        
%% Initialize variables 
P_lm=zeros(l_max+1,l_max+1,length(t));
u=sqrt(1-t.^2);
l=0:1:l_max;
%compute diagonal elements
P_lm(1,1,:)=1;
P_lm(2,2,:)=sqrt(3)*u;

for i=3:l_max+1
    P_lm(i,i,:)=sqrt((2*l(i)+1)/(2*l(i)))*u.*squeeze(P_lm(i-1,i-1,:));
end

%compute the other elements
for m=0:l_max-1 %for each degree
    P_lm(m+2,m+1,:)=sqrt(2*m+3)*t.*squeeze(P_lm(m+1,m+1,:));
    for i=m+2:l_max
        %i is l
        %m is m
        a=sqrt((2*i-1)*(2*i+1)/((i-m)*(i+m)));
        b=sqrt((2*i+1)*(i+m-1)*(i-m-1)/((2*i-3)*(i+m)*(i-m)));
        P_lm(i+1,m+1,:)=a*t.*squeeze(P_lm(i,m+1,:))-b*squeeze(P_lm(i-1,m+1,:));
    end
end
P_lm=P_lm;
end
