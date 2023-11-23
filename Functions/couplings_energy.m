%% Function used to compute the couplings for energy dissipation 
%% INPUT 
       % n_sol: degrees of modes involved in the solution \alpha
       % m_sol: orders involved in the solution \alpha
       % Nmax: maxium degree for which energy tidal patterns are computed
%% OUTPUT
       % EC: Matrix containing the couplings 
        % EC(n_mode1,n_mode2,n_mode3,n2_mode1, n2_mode2)
            % with n_mode1: mode 1 of the solution with n_s(n_mode1), m_s(n_mode1)
            % with n_mode2: mode 2 of the solution with n_s(n_mode1), m_s(n_mode1)
            % with n_mode3: mode of the energy spectra with n_en(n_mode1), m_en(n_mode1)
            % with n2_mode1: n2 of mode 1
            % with n2_mode2: n2 of mode 2
         
        % n_en: energy spectrum 

            

%% FUNCTION STARTS HERE
function [EC,n_en,m_en] = couplings_energy(n_sol,m_sol,Nmax)
%% compute n_en 
k=1;
for i=0:Nmax
    for j=-i:i
        n_en(k)=i; 
        m_en(k)=j; 
        k=k+1; 
    end   
end
%% initialize matrix
tic
EC=zeros(length(n_sol),length(n_sol),length(n_en),6,6);
for i1=1:size(EC,1) %loop modes 1
    for i2=1:size(EC,2)
        %(1) loop over modes that result from coupling
        % Selection rule (2). Triangular inequality.
        n_aux=[abs(n_sol(i1)-n_sol(i2)) n_sol(i1)+n_sol(i2)]; % Changed to allow for lower Nmax
        % Seclction rule (7)
        m_aux=m_sol(i1)+m_sol(i2); % 
        ind_resulting_modes=find(n_en>=n_aux(1) & n_en<=n_aux(end) & m_en==m_aux);        
        %(2) loop over modes that result from coupling 
        for i3=1:length(ind_resulting_modes)
            na=n_sol(i1);
            ma=m_sol(i1);
            nb=n_sol(i2);
            mb=m_sol(i2); % 
            nc=n_en(ind_resulting_modes(i3));
            mc=m_en(ind_resulting_modes(i3));
            %(3) loop over n2 coefficients
            la=[2 2 2 2 2 0];
            lb=[2 2 2 2 2 0];
            na2=[na-2:na+2 na];
            nb2=[nb-2:nb+2 nb];
            for i4=1:length(na2)
                for i5=1:length(nb2)
                    %compute couplings here
                    % Selection rule (4)
                    % Selection rule (3)
                    if mod(na2(i4)+nb2(i5)+nc,2)==0 && nb2(i5)>=abs(na2(i4)-nc) && nb2(i5)<=abs(na2(i4)+nc) 
                       EC(i1,i2,ind_resulting_modes(i3),i4,i5)=couplings_coefficient(na,na2(i4),la(i4),ma,nb,nb2(i5),lb(i5),mb,nc,mc); 
                    end
                    k=k+1; 
                end
            end
        end
    end
end
toc
%% find the non-zero coefficients 
k=1;
non_zero_index=[]; 
n_en2=[];
m_en2=[];
for i=1:length(n_en)
    EC_deg=abs(EC(:,:,i,:,:));
    max_c=max(EC_deg(:));
    if max_c>0
        n_en2(k)=n_en(i);
        m_en2(k)=m_en(i); 
        non_zero_index(k)=i;
        k=k+1; 
    end
end
EC=EC(:,:,non_zero_index,:,:); 
n_en=n_en2; 
m_en=m_en2;
%% make plot of solution and energy spectra degrees and order?
% k=1;
% for i=0:Nmax
%     for j=-i:i
%         n_v(k)=i; 
%         m_v(k)=j; 
%         k=k+1; 
%     end   
% end
% Acoup=zeros(length(n_v),2);
% for i=1:length(n_v)
%     in1=find(n_v(i)==n_sol & m_v(i)==m_sol);
%     in2=find(n_v(i)==n_en & m_v(i)==m_en);
%     if isempty(in1)==0
%         Acoup(i,1)=1;
%     end
%     if isempty(in2)==0
%         Acoup(i,2)=1;
%     end
% end
% Acoup(:,end+1)=0;
% Acoup(end+1,:)=0;
% hh=figure
% set(hh, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.3, 0.9]);
% pcolor(Acoup)
% colormap(flip(cmocean('gray')))
% hold on
% i=1; 
%     j_l=[];
%     i_L=[];
%     k=1; 
%     for n=0:Nmax
%         for m=-n:1:n
%           labels{i}=['(' num2str(n) ',' num2str(m) ')'];
%           m_v(i)=m;
%           n_v(i)=n;
%           if m==0 
%               j_l=[j_l i-0.5+1];
%               labels_P{k}=labels{i};
%               k=k+1; 
%           end
%           if m==n
%               i_L=[i_L i+1]; 
%           end
%           i=i+1;   
%         end
%     end
% % major grid lines, separating degree
% for k=1:length(i_L)
%     plot([i_L(k) i_L(k)],[0 i-1],'LineWidth',2,'color','k');
%     hold on
%     plot([0 i-1],[i_L(k) i_L(k)],'LineWidth',2,'color','k');
%     hold on
% end
% % minor grid lines, separating degree
% for k=1:i-1
%     plot([k k],[0 i-1],'LineWidth',0.1,'color',[0.5 0.5 0.5]);
%     hold on
%     plot([0 i-1],[k k],'LineWidth',0.1,'color',[0.5 0.5 0.5]);
%     hold on
% end
% %grid for the coefficients
% plot([i+3 i+3],[0 i-1],'LineWidth',2,'color','k');
% hold on 
% plot([2 2],[0 length(n_v)+1],'LineWidth',2,'color','k');
% hold on 
% plot([3 3],[0 length(n_v)+1],'LineWidth',2,'color','k');
% hold on
% plot([i+2 i+2],[0 i],'LineWidth',2,'color','k');
% for k=1:length(i_L)
%     plot([i+2 i+3],[i_L(k) i_L(k)],'LineWidth',2,'color','k');
%     hold on
% end
% box on
% xlim([1 3]);
% yticks(j_l)
% yticklabels(labels_P)
% xticks([1.5,2.5])
% xticklabels({'Response','Energy'})
% % axesH = axes;
% % axesH.XAxis.MinorTick       = 'on';
% % axesH.XAxis.MinorTickValues =0:1:i;
% ylim([1 i]);
% %xlabel('Equation')
% ylabel('Coefficients')
end

%% explicit computation of the coupling coefficients
function [C] = couplings_coefficient(na,na2,la,ma,nb,nb2,lb,mb,nc,mc)
% la=2;
% lb=2;
na1v=[na-1 na na+1];
nb1v=[nb-1 nb nb+1];
C=0;
Ctest=0;
for i=1:length(na1v)
    na1=na1v(i); 
    for j=1:length(nb1v)
        nb1=nb1v(j);
        Lama=sqrt((2*la+1)*(2*na1+1));
        Lamb=sqrt((2*lb+1)*(2*nb1+1));
        CC=(-1)^(mc+nb+nb2)*...
        sqrt((2*na2+1)*(2*na1+1)*(2*na+1))*...
        sqrt((2*nb2+1)*(2*nb1+1)*(2*nb+1))*...
        sqrt(2*nc+1)*...
        Wigner3j(na2, nb2, nc, 0, 0, 0)*...
        Wigner3j(na,nb,nc,ma,mb,-mc)*...
        Wigner9j(na,na1,1,nc,na2,nb2,nb,1,nb1);
        Caux=(-1)^(na+na2+la+nb+nb2+lb)*Lama*Lamb*Wigner6j(1, la, 1, na, na1, na2)*Wigner6j(1, lb, 1, nb, nb1, nb2)*CC;
        C=C+Caux; 
    % test for degree 0 coupling
%         if nc==0 && na==nb && na2==nb2
%             CC=(-1)^(na1+nb1)*sqrt((2*na1+1))*sqrt((2*nb2+1)*(2*nb1+1)*(2*nb+1))*...
%             Wigner3j(na2, nb2, nc, 0, 0, 0)*...
%             Wigner3j(na,nb,nc,ma,mb,0)*...
%             Wigner6j(na1,nb,1,nb1,nb2,1);
%             Caux=(-1)^(na+na2+la+nb+nb2+lb)*Lama*Lamb*Wigner6j(1, la, 1, na, na1, na2)*Wigner6j(1, lb, 1, nb, nb1, nb2)*CC;
%             Ctest=Ctest+Caux;
%         end
    end
end

% if nc==0 && na==nb && na2==nb2
%     C-Ctest
% end
end

