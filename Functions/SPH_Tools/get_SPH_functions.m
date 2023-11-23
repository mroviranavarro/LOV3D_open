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
