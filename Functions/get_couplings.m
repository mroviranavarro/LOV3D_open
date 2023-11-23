%% GET_COUPLINGS
% AUTHOR: M. Rovira-Navarro & A. Veenstra 
% USE: function used to get the coupling coefficients
%% INPUT 
    % max_order: maximum perturbation order, used to obtain the modes of the solution 
    % variations: contians the degree, variations(:,1), and order, variations(:,2) of the lateral variations. The length is Nreo. 
    % Forcings: spectra of the foricng. 
        % Forcing.Td: forcing period
        % Forcing.n: degree of the forcing 
        % Forcing.m: order of the forcing 
        % Forcing.F: amplitude of the component 
%% OUTPUT
% Couplings 
    % Couplings.n_s: degree of modes that intervene in the solution
    % Couplings.m_s: order of modes that intervene in the solution
    % Couplings.order: order of the coupling 
    % Couplings.Coup: matrix containing coupling coefficients 
        % Coup(Nsol,Nsol,27,Nreo)
            % Coup(i,j,1,k)  (T^ma_{na,na,0}Y_nb^mb,T^mc_{nc,nc,0})        
            % Coup(i,j,2,k)  (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc,2})
            % Coup(i,j,3,k)  (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc,2})
            % Coup(i,j,4,k)  (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc,2})
            % Coup(i,j,5,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc,2})
            % Coup(i,j,6,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc,2})        
            % Coup(i,j,7,k)   (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc-2,2})
            % Coup(i,j,8,k)   (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc-2,2})
            % Coup(i,j,9,k)   (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc-2,2})
            % Coup(i,j,10,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc-2,2})
            % Coup(i,j,11,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc-2,2})        
            % Coup(i,j,12,k)  (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc+2,2})
            % Coup(i,j,13,k)  (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc+2,2})
            % Coup(i,j,14,k)  (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc+2,2})
            % Coup(i,j,15,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc+2,2})
            % Coup(i,j,16,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc+2,2})        
            % Coup(i,j,17,k)  (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc-1,2})
            % Coup(i,j,18,k)  (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc-1,2})
            % Coup(i,j,19,k)  (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc-1,2})
            % Coup(i,j,20,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc-1,2})
            % Coup(i,j,21,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc-1,2})        
            % Coup(i,j,22,k)  (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,23,k)  (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,24,k)  (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,25,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,26,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,27,k) 0 if all coupling coefficents are 0

%% FUNCTION STARTS HERE 
function [Couplings]=get_couplings(max_order,variations,Forcing,varargin) 
%% OPTIONAL INPUTS
verbose=0;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose=1;
        varargin{k}=[];
    end
end
%% GET MODES THAT ARE COUPLED USING SELECTION RULES
% loop over rheology
modes=[];
rheo=unique((variations(:,[1 2])),'rows');
Nreo=size(variations,1);

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

Couplings.n_s=modes4(:,1);
Couplings.m_s=modes4(:,2);
Couplings.order=modes4(:,3);
Nsol=length(Couplings.n_s);
%% GET COUPLING COEFFICIENTS 
if verbose==1
     disp("Getting in the loop to compute coupling coefficients...")
     disp([ num2str(Nsol) ' Modes'])
     disp([ num2str(Nreo) ' Rheology'])
      disp([ num2str((Nreo*Nsol^2)) ' Couplings to be computed'])
     tic
end
Cp=zeros(Nsol,Nsol,27,Nreo);
jj=1;
for ireo=1:Nreo
    nb=variations(ireo,1);
    mb=variations(ireo,2);
    for imod1=1:Nsol
        n=Couplings.n_s(imod1);
        m=Couplings.m_s(imod1);
        for imod2=1:Nsol
            na=Couplings.n_s(imod2);
            ma=Couplings.m_s(imod2);
            % check if the coefficient is non-zero
            % abs(n-nb):n+nb 
            if na>=abs(n-nb) && na<=abs(n+nb) && ma==m-mb
                Cp(imod1,imod2,:,ireo)=coupling_coefficients(n,m,na,ma,nb,mb);              
            end
            jj=jj+1;
            if verbose==1
                clc
                disp([ num2str(Nsol) ' Modes'])
                disp([ num2str(Nreo) ' Rheology'])
                disp([ num2str((Nreo*Nsol^2)) ' Couplings to be computed'])
                disp([ num2str(jj/(Nreo*Nsol^2)*100) ' % Couplings computed'])
                toc
            end
        end
    end
end
Couplings.Coup=Cp;
end
%% FUNCTION USED TO GET NEXT SET OF COUPLINGS
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
%% FUNCTION TO GET THE COUPLINGS 
function Cp=coupling_coefficients(n,m,na,ma,nb,mb)
Cp(1,1,1)=couplingT(n,n,0,m,na,na,0,ma,nb,mb);                  
%l=2 n2=n; 
Cp(1,1,2)=couplingT(n,n,2,m,na,na,2,ma,nb,mb);
Cp(1,1,3)=couplingT(n,n,2,m,na,na-2,2,ma,nb,mb);
Cp(1,1,4)=couplingT(n,n,2,m,na,na+2,2,ma,nb,mb);
Cp(1,1,5)=couplingT(n,n,2,m,na,na-1,2,ma,nb,mb);
Cp(1,1,6)=couplingT(n,n,2,m,na,na+1,2,ma,nb,mb);                    
%l=2 n2=n-2; 
Cp(1,1,7)=couplingT(n,n-2,2,m,na,na,2,ma,nb,mb);
Cp(1,1,8)=couplingT(n,n-2,2,m,na,na-2,2,ma,nb,mb);
Cp(1,1,9)=couplingT(n,n-2,2,m,na,na+2,2,ma,nb,mb);
Cp(1,1,10)=couplingT(n,n-2,2,m,na,na-1,2,ma,nb,mb);
Cp(1,1,11)=couplingT(n,n-2,2,m,na,na+1,2,ma,nb,mb);                    
%l=2 n2=n+2; 
Cp(1,1,12)=couplingT(n,n+2,2,m,na,na,2,ma,nb,mb);
Cp(1,1,13)=couplingT(n,n+2,2,m,na,na-2,2,ma,nb,mb);
Cp(1,1,14)=couplingT(n,n+2,2,m,na,na+2,2,ma,nb,mb);
Cp(1,1,15)=couplingT(n,n+2,2,m,na,na-1,2,ma,nb,mb);
Cp(1,1,16)=couplingT(n,n+2,2,m,na,na+1,2,ma,nb,mb);                    
%l=2 n2=n-1; 
Cp(1,1,17)=couplingT(n,n-1,2,m,na,na,2,ma,nb,mb);
Cp(1,1,18)=couplingT(n,n-1,2,m,na,na-2,2,ma,nb,mb);
Cp(1,1,19)=couplingT(n,n-1,2,m,na,na+2,2,ma,nb,mb);
Cp(1,1,20)=couplingT(n,n-1,2,m,na,na-1,2,ma,nb,mb);
Cp(1,1,21)=couplingT(n,n-1,2,m,na,na+1,2,ma,nb,mb);                    
%l=2 n2=n+1; 
Cp(1,1,22)=couplingT(n,n+1,2,m,na,na,2,ma,nb,mb);
Cp(1,1,23)=couplingT(n,n+1,2,m,na,na-2,2,ma,nb,mb);
Cp(1,1,24)=couplingT(n,n+1,2,m,na,na+2,2,ma,nb,mb);
Cp(1,1,25)=couplingT(n,n+1,2,m,na,na-1,2,ma,nb,mb);
Cp(1,1,26)=couplingT(n,n+1,2,m,na,na+1,2,ma,nb,mb);                                         
if max(abs(squeeze(Cp(:))))>0
    Cp(1,1,27)=1;
else
    Cp(1,1,27)=0;
end
end

function C=couplingT(nc,nc2,lc,mc,na,na2,la,ma,nb,mb)
na1v=[na-1 na na+1];
nc1v=[nc-1 nc nc+1];
C=0; 
for i=1:length(na1v)
        na1=na1v(i);
    for j=1:length(nc1v)
        nc1=nc1v(j);
        Lama=sqrt((2*la+1)*(2*na1+1));
        Lamc=sqrt((2*lc+1)*(2*nc1+1));
        Caux=(-1)^(na+na2+la+nc+nc2+lc)*Lama*Lamc*Wigner6j(1, la, 1, na, na1, na2)*Wigner6j(1, lc, 1, nc, nc1, nc2)*couplingY(nc,mc,nc1,nc2,na,ma,na1,na2,nb,mb);
        C=C+Caux; 
    end
end
end

function [C] = couplingY(nc,mc,nc1,nc2,na,ma,na1,na2,nb,mb)
Lam=sqrt((2*nc2+1)*(2*nc1+1)*(2*nc+1)*...
        (2*na2+1)*(2*na1+1)*(2*na+1)*...
        (2*nb+1));   
C=(-1)^(na+na1+nc+nc1+mc)*Lam*...
       Wigner6j(na, na1, 1, nc1, nc, nb)*...
       Wigner6j(na1, na2, 1, nc2, nc1, nb)*...
       Wigner3j(na2, nc2, nb, 0, 0, 0)*...
       Wigner3j(na,nc,nb,ma,-mc,mb);
end