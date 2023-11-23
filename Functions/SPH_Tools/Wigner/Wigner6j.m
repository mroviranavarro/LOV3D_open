function wig = Wigner6j(j1,j2,j3,j4,j5,j6,ifs,ifcb)
% Computes the Wigner 6j symbols
%    wig = {j1 j2 j3
%           j4 j5 j6}
% or the recoupling matrix element
%    wig = <(j1,j2)j3,j4,j5|j1,(j2,j4)j6,j5>
% using the Racah formula.
%
%     USAGE
% wig = Wigner6j(j1,j2,j3,j4,j5,j6,ifs,ifcb);
% Wigner6j(j1,j2,j3,j4,j5,j6)
%
%     INPUT
%     mandatory arguments:
% j1, j2, j3, j4, j5, j6
%      - the sets of the angular momenta; arrays of the same sizes;
%     optional arguments:
% ifs  - the switch of the computational methods;
%        although the central part of all the computations is based on the Racah formula, some details differ:
%      = 0 - the symbolic (presumably accurate) computation with the double-precision output;
%      =-1 - the same symbolic computation as ifs=0 but with the symbolic output (accurate square root/rational-type equation, simplified);
%      =-2 - the same as ifs=-1 but without a final simplification of the symbolic expression, which can be time-consuming sometimesl
%      = 1, =2 (default, recommended), =3 - numeric double-precision algorithms (see REMARKS below);
%      all other input values of ifs are set to the closest from the above list.
% ifcb - if exists and is true, switches to computing the coupling matrix elements instead of the 6j-symbols;
%        default ifcb=false (6j-symbols are computed);
%     OUTPUT
% wig  - the resulting values of the 6j-symbols (ifcb=false) or the coupling matrix element (ifcb=true) in either numeric double-precision or symbolic form (see the input parameter ifs);
%        array of the same size as j1.
%     REMARKS:
%     most of the estimates below is based on extensive numerical tests within a range of the quantum nubers <=1000
% (a) the numeric algorithms (ifs>0) are usually much faster than the symbolic algorithm (ifs<=0);
% (b) the accuracy of the symbolic algorithm remains the uttermost possible;
% (c) the accuracy of the numeric algorithms can worsen for big quantum numbers;
% (d) all the "numeric" algorithms switch automatically to the symbolic computations
%     as soon as the numeric overflow occurs, thereby improving the accuracy
%     but slowing down the computations for some big quantum numbers;
% (e) in cases with at least one of the j quantum numbers <=2, the explicit accurate equations
%     are applied in all the algorithms ensuring the highest possible accuracy;
% (f) for relatively small quantum numbers (up to ~20) all the algorithms provide
%     approximately equivalent results with a reasonably high accuracy;
% (g) the algorithms ifs=1 and ifs=2 ensure a reasonably good accuracy; the worst registered cases
%     for both algorithms (occuring before switching to the symbolic regime) was;
%     Wigner6j(41.5,52,43.5,36,38.5,46)=-0.00029 with the biggest relative inaccuracy of 3.6e-10
% (h) the algorithm ifs=2 (comparing to ifs=1) proceeds in a wider range of quantum numbers numerically
%     before switching to the symbolic regime (i.e., can work faster with big quantum numbers);
%     in our experimental probing it proved to ensure a good enough accuracy and was chosen as a default one;
% (i) the algorithm ifs=3 (comparing to ifs=2) proceeds in a wider range of quantum numbers numerically
%     before switching to the symbolic regime (i.e., can work faster with big quantum numbers)
%     but for some combinations of big quantum numbers, completely wrong results are observed:
%     e.g., Wigner6j(227,230.5,210.5,249.5,235,232,3);
%     thereby, we do not recommend using the algorithm ifs=3 at all, and only keep it here in view of the completeness.
%
%     LINKS (DEFINITIONS)
% (1) https://en.wikipedia.org/wiki/6-j_symbol
% (2) R. N. Zare, Angular Momentum: Understanding Spatial Aspects in Chemistry and
%     Physics, John Wiley & Sons, New York砲hichester烹risbane傍oronto亡igapore, 1988.
%     URL: http://eu.wiley.com/WileyCDA/WileyTitle/productCd-0471858927.html
%
% VERSION of January, 2020
%
% Programmer: V. B. Sovkov
% St. Petersburg State University
% Shanxi University
%%

try


    NJ = min([numel(j1) numel(j2) numel(j3) numel(j4) numel(j5) numel(j6)]);
    Nj1S=size(j1);
    j1 = double(reshape(j1(1:NJ),NJ,1));
    j2 = double(reshape(j2(1:NJ),NJ,1));
    j3 = double(reshape(j3(1:NJ),NJ,1));
    j4 = double(reshape(j4(1:NJ),NJ,1));
    j5 = double(reshape(j5(1:NJ),NJ,1));
    j6 = double(reshape(j6(1:NJ),NJ,1));

    if nargin<7 || isempty(ifs) || ~isnumeric(ifs) && ~islogical(ifs) || ifs(1)>0 && ifs(1)<=1.5 % choose the algorithm
        ifs = 2; % default
    elseif ifs(1)>0
        ifs = min(round(ifs(1)),3);
    else
        ifs=max(round(ifs(1)),-2);
    end

    NK = if3jc(j1,j2,j3) & if3jc(j1,j5,j6) & if3jc(j3,j4,j5) & if3jc(j2,j4,j6);
    if all(~NK)
        if ifs>=0
            wig=zeros(NJ,1);
        else
            wig=sym(zeros(NJ,1));
        end
        if prod(Nj1S)==NJ
            reshape(wig,Nj1S);
        end
        return;
    end

    if ifs>0
        wig=zeros(NJ,1);
    else
        wig=sym(zeros(NJ,1));
    end

    ifcb = nargin>7 && ~isempty(ifcb) && ifcb(1); % if the Clebsch-Gordan coefficient instead of the 3j-symbols are needed?
    if ifcb
        J1=j1;
        J2=j2;
        J3=j3;
        J4=j4;
        J5=j5;
        J6=j6;
    end

    for k=1:NJ
        jup=[j1(k),j2(k),j3(k)]';
        jdw=[j4(k),j5(k),j6(k)]';
        if min(jdw)<min(jup)
            [~,jup,jdw]=ArranA(jdw,jup);
            jj=jdw(3);
            jdw(3)=jup(3);
            jup(3)=jj;
        end
        [~,jup,jdw]=ArranA(jup,jdw);
        while jup(2)>jdw(2)
            jj=jdw(2:3);
            jdw(2:3)=jup(2:3);
            jup(2:3)=jj;
            [~,jup,jdw]=ArranA(jup,jdw);
        end
        j1(k)=jup(1);
        j2(k)=jup(2);
        j3(k)=jup(3);
        j4(k)=jdw(1);
        j5(k)=jdw(2);
        j6(k)=jdw(3);
    end
    if NJ>1
        [ind0,j1,j2,j3,j4,j5,j6] = ...
            ArranA ( reshape(j1(1:NJ),NJ,1),reshape(j2(1:NJ),NJ,1),reshape(j3(1:NJ),NJ,1),reshape(j4(1:NJ),NJ,1),reshape(j5(1:NJ),NJ,1),reshape(j6(1:NJ),NJ,1));
        NK = NK(ind0);

        kk = find ( j1(1:NJ-1)==j1(2:NJ) & j2(1:NJ-1)==j2(2:NJ) & j3(1:NJ-1)==j3(2:NJ) & j4(1:NJ-1)==j4(2:NJ) & j5(1:NJ-1)==j5(2:NJ) & j6(1:NJ-1)==j6(2:NJ) );
        if ~isempty(kk)
            NK(kk)=0;
        end
    else
        kk=[];
    end

    NKs=[]; % for pointing to erronious numeric results in order to switch them to the symbolic computations

%%
%   special cases - faster and more accurate computation than the general case below

    if any(NK & j1<=2) % the range of the special cases currently included

%   case j = 0
        kp = find(NK & j1==0);
        if ~isempty(kp)
            if ifs>0
                wig(kp) = (-1).^(j2(kp)+j4(kp)+j5(kp)) ./ sqrt((2*j2(kp)+1).*(2*j5(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^(j2(kp)+j4(kp)+j5(kp)) ./ sqrt(sym((2*j2(kp)+1).*(2*j5(kp)+1)));
            end
            NK(kp)=0;
        end
%   case   j = 1/2
        kp = find(NK & j1==1/2 & j5<j6);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt((s-2*j5(kp))./(2*j5(kp)+1)./(2*j5(kp)+2).*(s-2*j3(kp)+1)./(2*j3(kp))./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym((s-2*j5(kp))./(2*j5(kp)+1)./(2*j5(kp)+2).*(s-2*j3(kp)+1)./(2*j3(kp))./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1/2); %  & j5>j6
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt((s+1)./(2*j5(kp))./(2*j5(kp)+1).*(s-2*j4(kp))./(2*j3(kp))./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym((s+1)./(2*j5(kp))./(2*j5(kp)+1).*(s-2*j4(kp))./(2*j3(kp))./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
%   case   j = 1
        kp = find(NK & j1==1 & j5>j6 & j3>j2); % & j3>j2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(s./(2*j5(kp)-1)./(2*j5(kp)).*(s+1)./(2*j5(kp)+1).*(s-2*j4(kp)-1)./(2*j3(kp)-1).*(s-2*j4(kp))./(2*j3(kp))./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(s./(2*j5(kp)-1)./(2*j5(kp)).*(s+1)./(2*j5(kp)+1).*(s-2*j4(kp)-1)./(2*j3(kp)-1).*(s-2*j4(kp))./(2*j3(kp))./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j5<j6 & j3>j2); % & j3>j2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt((s-2*j5(kp)-1)./(2*j5(kp)+1)./(2*j5(kp)+2).*(s-2*j5(kp))./(2*j5(kp)+3)./(2*j3(kp)-1).*(s-2*j3(kp)+1)./(2*j3(kp)).*(s-2*j3(kp)+2)./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym((s-2*j5(kp)-1)./(2*j5(kp)+1)./(2*j5(kp)+2).*(s-2*j5(kp))./(2*j5(kp)+3)./(2*j3(kp)-1).*(s-2*j3(kp)+1)./(2*j3(kp)).*(s-2*j3(kp)+2)./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j5==j6 & j3>j2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(2*(s+1)./(2*j5(kp))./(2*j5(kp)+1).*(s-2*j4(kp))./(2*j5(kp)+2)./(2*j3(kp)-1).*(s-2*j5(kp))./(2*j3(kp)).*(s-2*j3(kp)+1)./(2*j3(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(2*(s+1)./(2*j5(kp))./(2*j5(kp)+1).*(s-2*j4(kp))./(2*j5(kp)+2)./(2*j3(kp)-1).*(s-2*j5(kp))./(2*j3(kp)).*(s-2*j3(kp)+1)./(2*j3(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j5<j6 & j3==j2);
        if ~isempty(kp)
            s=j2(kp)+j4(kp)+j6(kp);
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(2*(s+1)./(2*j2(kp))./(2*j2(kp)+1).*(s-2*j4(kp))./(2*j2(kp)+2)./(2*j6(kp)-1).*(s-2*j2(kp))./(2*j6(kp)).*(s-2*j6(kp)+1)./(2*j6(kp)+1));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(2*(s+1)./(2*j2(kp))./(2*j2(kp)+1).*(s-2*j4(kp))./(2*j2(kp)+2)./(2*j6(kp)-1).*(s-2*j2(kp))./(2*j6(kp)).*(s-2*j6(kp)+1)./(2*j6(kp)+1)));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j5==j6  & j2==j3); % & j2==j3
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
            if ifs>0
                wig(kp) = 2*(-1).^s .* (j4(kp).*(j4(kp)+1)-j5(kp).*(j5(kp)+1)-j3(kp).*(j3(kp)+1)) ./ sqrt((2*j5(kp)).*(2*j5(kp)+1).*(2*j5(kp)+2).*(2*j3(kp)).*(2*j3(kp)+1).*(2*j3(kp)+2));
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .* sym(j4(kp).*(j4(kp)+1)-j5(kp).*(j5(kp)+1)-j3(kp).*(j3(kp)+1)) ./ sqrt(sym((2*j5(kp)).*(2*j5(kp)+1).*(2*j5(kp)+2).*(2*j3(kp)).*(2*j3(kp)+1).*(2*j3(kp)+2)));
            end
            NK(kp)=0;
        end
%   case   j = 3/2
        kp = find(NK & j1==3/2 & j3-j2==3/2 & j5-j6==3/2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+1;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       (s-1).*s.*(s+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)-2).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       (s-1).*s.*(s+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)-2).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==3/2 & j5-j6==1/2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+2;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       3*s.*(s+1).*(s-2*j4(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       3*s.*(s+1).*(s-2*j4(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kp = find(NK & j1==3/2 & j3-j2==1/2 & j5-j6==3/2);
        if ~isempty(kp)
            s=j4(kp)+j3(kp)+j5(kp);
			A1=2*j3(kp)+2;
			A2=2*j5(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       3*s.*(s+1).*(s-2*j4(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)).*(s-2*j3(kp)).*(s-2*j5(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       3*s.*(s+1).*(s-2*j4(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j4(kp)).*(s-2*j3(kp)).*(s-2*j5(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==1/2 & j5-j6==-3/2);
        if ~isempty(kp)
            s=j4(kp)+j2(kp)+j6(kp);
			A1=2*j2(kp)+3;
			A2=2*j6(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       3*(s+1).*(s-2*j4(kp)).*(s-2*j2(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       3*(s+1).*(s-2*j4(kp)).*(s-2*j2(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kp = find(NK & j1==3/2 & j3-j2==3/2 & j5-j6==-1/2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       3*(s+1).*(s-2*j4(kp)).*(s-2*j5(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       3*(s+1).*(s-2*j4(kp)).*(s-2*j5(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==3/2 & j5-j6==-3/2); % & j5-j6==-3/2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+4;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .* sqrt(...
                       (s-2*j5(kp)-2).*(s-2*j5(kp)-1).*(s-2*j5(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .* sqrt(sym(...
                       (s-2*j5(kp)-2).*(s-2*j5(kp)-1).*(s-2*j5(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    .* (s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j5-j6==1/2 & j3-j2==1/2); % & j3-j2==1/2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+2;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      (s+1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				    ) .* ...
                       (2*(s-2*j5(kp)).*(s-2*j3(kp))-(s+2).*(s-2*j4(kp)-1))...
				;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      (s+1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				    )) .* ...
                       sym(2*(s-2*j5(kp)).*(s-2*j3(kp))-(s+2).*(s-2*j4(kp)-1))...
				;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==1/2 & j5-j6==-1/2); % & j3-j2==1/2 & j5-j6==-1/2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      (s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				    ) .* ...
                       ((s-2*j5(kp)-1).*(s-2*j3(kp))-2*(s+2).*(s-2*j4(kp)))...
				;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      (s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)...
				    )) .* ...
                       sym((s-2*j5(kp)-1).*(s-2*j3(kp))-2*(s+2).*(s-2*j4(kp)))...
				;
            end
            NK(kp)=0;
        end
%   case   j = 2, j3-j2=2
        kp = find(NK & j1==2 & j3-j2==2 & j5-j6==2);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+1;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      (s-2).*(s-1).*s.*(s+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-3).*(s-2*j4(kp)-2).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      (s-2).*(s-1).*s.*(s+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-3).*(s-2*j4(kp)-2).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==2 & j5-j6==1);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+2;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      (s-1).*s.*(s+1).*(s-2*j4(kp)-2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-1).*(s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      (s-1).*s.*(s+1).*(s-2*j4(kp)-2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-1).*(s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==2 & j5==j6);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      6*s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j5(kp)-1).*(s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      6*s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j5(kp)-1).*(s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==2 & j5-j6==-1);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+4;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      (s+1).*(s-2*j4(kp)).*(s-2*j5(kp)-2).*(s-2*j5(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      (s+1).*(s-2*j4(kp)).*(s-2*j5(kp)-2).*(s-2*j5(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==2 & j5-j6==-2); % & j5-j6==-2
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+5;
			A2=2*j3(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      (s-2*j5(kp)-3).*(s-2*j5(kp)-2).*(s-2*j5(kp)-1).*(s-2*j5(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3).*(s-2*j3(kp)+4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      (s-2*j5(kp)-3).*(s-2*j5(kp)-2).*(s-2*j5(kp)-1).*(s-2*j5(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j3(kp)+1).*(s-2*j3(kp)+2).*(s-2*j3(kp)+3).*(s-2*j3(kp)+4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
%   case   j = 2, j3-j2=1
        kp = find(NK & j1==2 & j3-j2==1 & j5-j6==2);
        if ~isempty(kp)
            s=j4(kp)+j3(kp)+j5(kp);
			A1=2*j3(kp)+2;
			A2=2*j5(kp)+1;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      (s-1).*s.*(s+1).*(s-2*j4(kp)-2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-1).*(s-2*j4(kp)).*(s-2*j3(kp)).*(s-2*j5(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      (s-1).*s.*(s+1).*(s-2*j4(kp)-2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j4(kp)-1).*(s-2*j4(kp)).*(s-2*j3(kp)).*(s-2*j5(kp)+1)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & j5-j6==-2);
        if ~isempty(kp)
            s=j4(kp)+j2(kp)+j6(kp);
			A1=2*j2(kp)+4;
			A2=2*j6(kp)+1;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      (s+1).*(s-2*j4(kp)).*(s-2*j2(kp)-2).*(s-2*j2(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2).*(s-2*j6(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      (s+1).*(s-2*j4(kp)).*(s-2*j2(kp)-2).*(s-2*j2(kp)-1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2).*(s-2*j6(kp)+3)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & j5-j6==1);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+2;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = 4*(-1).^s .*...
                    sqrt(...
                      s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ) .*...
                      ( (j4(kp)+j5(kp)).*(j4(kp)-j5(kp)+1)-(j3(kp)-1).*(j3(kp)-j5(kp)+1) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 4*(-1).^s .*...
                    sqrt(sym(...
                      s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    )) .*...
                      sym( (j4(kp)+j5(kp)).*(j4(kp)-j5(kp)+1)-(j3(kp)-1).*(j3(kp)-j5(kp)+1) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & j5==j6);
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      6*(s+1).*(s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ) .*...
                      ( (j4(kp)+j5(kp)+1).*(j4(kp)-j5(kp))-j3(kp).^2+1 )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      6*(s+1).*(s-2*j4(kp)).*(s-2*j5(kp)).*(s-2*j3(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    )) .*...
                      sym( (j4(kp)+j5(kp)+1).*(j4(kp)-j5(kp))-j3(kp).^2+1 )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & j5-j6==-1); % & j5-j6==-1
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+4;
			A2=2*j3(kp)+2;
            if ifs>0
                wig(kp) = 4*(-1).^s .*...
                    sqrt(...
                      (s-2*j5(kp)-1).*(s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ) .*...
                      ( (j4(kp)+j5(kp)+2).*(j4(kp)-j5(kp)-1)-(j3(kp)-1).*(j5(kp)+j3(kp)+2) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 4*(-1).^s .*...
                    sqrt(sym(...
                      (s-2*j5(kp)-1).*(s-2*j5(kp)).*(s-2*j3(kp)+1).*(s-2*j3(kp)+2)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    )) .*...
                      sym( (j4(kp)+j5(kp)+2).*(j4(kp)-j5(kp)-1)-(j3(kp)-1).*(j5(kp)+j3(kp)+2) )...
                    ;
            end
            NK(kp)=0;
        end
%   case   j = 2, j3==j2
        kp = find(NK & j1==2 & j3==j2 & j5-j6==-2);
        if ~isempty(kp)
            s=j4(kp)+j2(kp)+j6(kp);
			A1=2*j2(kp)+3;
			A2=2*j6(kp)+1;
            if ifs>0
                wig(kp) = (-1).^s .*...
                    sqrt(...
                      6*s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j2(kp)-1).*(s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = (-1).^s .*...
                    sqrt(sym(...
                      6*s.*(s+1).*(s-2*j4(kp)-1).*(s-2*j4(kp))...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
					.* (s-2*j2(kp)-1).*(s-2*j2(kp)).*(s-2*j6(kp)+1).*(s-2*j6(kp)+2)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ));
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3==j2 & j5-j6==-1);
        if ~isempty(kp)
            s=j4(kp)+j2(kp)+j6(kp);
			A1=2*j2(kp)+3;
			A2=2*j6(kp)+2;
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    sqrt(...
                      6*(s+1).*(s-2*j4(kp)).*(s-2*j2(kp)).*(s-2*j6(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    ) .*...
                      ( (j4(kp)+j2(kp)+1).*(j4(kp)-j2(kp))-j6(kp).^2+1 )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sqrt(sym(...
                      6*(s+1).*(s-2*j4(kp)).*(s-2*j2(kp)).*(s-2*j6(kp)+1)...
                    ./ A1./(A1-1)./(A1-2)./(A1-3)./(A1-4)...
                    ./ A2./(A2-1)./(A2-2)./(A2-3)./(A2-4)...
				    )) .*...
                      sym( (j4(kp)+j2(kp)+1).*(j4(kp)-j2(kp))-j6(kp).^2+1 )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3==j2 & j5==j6); % & j3==j2 & j5==j6
        if ~isempty(kp)
            s=j3(kp)+j4(kp)+j5(kp);
			A1=2*j5(kp)+3;
			A2=2*j3(kp)+3;
			C =j4(kp).*(j4(kp)+1)-j5(kp).*(j5(kp)+1)-j3(kp).*(j3(kp)+1);
            if ifs>0
                wig(kp) = 2*(-1).^s .*...
                    (3*C.*(C+1)-4*j5(kp).*(j5(kp)+1).*j3(kp).*(j3(kp)+1)) ./...
                    sqrt(...
                    A1.*(A1-1).*(A1-2).*(A1-3).*(A1-4)...
                    .* A2.*(A2-1).*(A2-2).*(A2-3).*(A2-4)...
				    );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = 2*(-1).^s .*...
                    sym(3*C.*(C+1)-4*j5(kp).*(j5(kp)+1).*j3(kp).*(j3(kp)+1)) ./...
                    sqrt(sym(...
                    A1.*(A1-1).*(A1-2).*(A1-3).*(A1-4)...
                    .* A2.*(A2-1).*(A2-2).*(A2-3).*(A2-4)...
				    ));
            end
            NK(kp)=0;
        end

    end
%%

    NK = find(NK);

    if ~isempty(NK)
        if ifs==3
            C =     tc(j1(NK),j2(NK),j3(NK),ifs)  + tc(j1(NK),j5(NK),j6(NK),ifs)  + tc(j4(NK),j2(NK),j6(NK),ifs)  + tc(j4(NK),j5(NK),j3(NK),ifs);
        else
            C =     tc(j1(NK),j2(NK),j3(NK),ifs) .* tc(j1(NK),j5(NK),j6(NK),ifs) .* tc(j4(NK),j2(NK),j6(NK),ifs) .* tc(j4(NK),j5(NK),j3(NK),ifs);
        end
        k0 = find(isinf(C) | isnan(C));
        if ~isempty(k0)
            NKs = [NKs;NK(k0)];
            NK(k0) = [];
            C(k0)=[];
        end
        
        if ~isempty(NK)

            t1 = max ( [j1(NK) + j2(NK) + j3(NK), j1(NK) + j5(NK) + j6(NK), j4(NK) + j2(NK) + j6(NK), j4(NK) + j5(NK) + j3(NK)] , [] , 2 );
            t2 = min ( [j1(NK) + j2(NK) + j4(NK) + j5(NK), j2(NK) + j3(NK) + j5(NK) + j6(NK), j3(NK) + j1(NK) + j6(NK) + j4(NK)] , [] , 2 );

            if ifs>0 % numeric
                if ifs==3
                    for k=numel(NK):-1:1 % t below can be of different lengths for the different k! Cannot apply the element-wise operations in this sense!
                        kN = NK(k);
                        t = (t1(k):t2(k))';
                        wig(kN) = ((-1).^t)' * exp ( gammaln([...
                            (t+2),(t-j1(kN)-j2(kN)-j3(kN)+1),(t-j1(kN)-j5(kN)-j6(kN)+1),(t-j4(kN)-j2(kN)-j6(kN)+1),(t-j4(kN)-j5(kN)-j3(kN)+1),(j1(kN)+j2(kN)+j4(kN)+j5(kN)-t+1),(j2(kN)+j3(kN)+j5(kN)+j6(kN)-t+1),(j3(kN)+j1(kN)+j6(kN)+j4(kN)-t+1)...
                            ]) * [1;-1;-1;-1;-1;-1;-1;-1] + C(k) );
                    end
                else % ifs==1 || ifs==2
                    for k=numel(NK):-1:1 % t below can be of different lengths for the different k! Cannot apply the element-wise operations in this sense!
                        kN = NK(k);
                        t = (t1(k):t2(k))';
                        wig(kN) = (-1).^t' * (...
                            factorial(t+1)...
                            ./factorial(t-j1(kN)-j2(kN)-j3(kN))...
                            * C(k)...
                            ./factorial(t-j1(kN)-j5(kN)-j6(kN))...
                            ./factorial(t-j4(kN)-j2(kN)-j6(kN))...
                            ./factorial(t-j4(kN)-j5(kN)-j3(kN))...
                            ./factorial(j1(kN)+j2(kN)+j4(kN)+j5(kN)-t)...
                            ./factorial(j2(kN)+j3(kN)+j5(kN)+j6(kN)-t)...
                            ./factorial(j3(kN)+j1(kN)+j6(kN)+j4(kN)-t)...
                            );
                    end
                end
                k0 = find(isinf(wig(NK)) | isnan(wig(NK)));
                if ~isempty(k0)
                    NKs = [NKs;NK(k0)];
                end
            else % symbolic
                for k=numel(NK):-1:1 % t below can be of different lengths for the different k! Cannot apply the element-wise operations in this sense!
                    kN = NK(k);
                    t = (t1(k):t2(k))';
                    wig(kN) = (-1).^t' * ( ...
                        factorial(sym(t+1))...
                        ./factorial(sym(t-j1(kN)-j2(kN)-j3(kN)))...
                        * C(k)...
                        ./factorial(sym(t-j1(kN)-j5(kN)-j6(kN)))...
                        ./factorial(sym(t-j4(kN)-j2(kN)-j6(kN)))...
                        ./factorial(sym(t-j4(kN)-j5(kN)-j3(kN)))...
                        ./factorial(sym(j1(kN)+j2(kN)+j4(kN)+j5(kN)-t))...
                        ./factorial(sym(j2(kN)+j3(kN)+j5(kN)+j6(kN)-t))...
                        ./factorial(sym(j3(kN)+j1(kN)+j6(kN)+j4(kN)-t))...
                        );
                end
            end
        end
    end
%%
%   post-processing

    if~ifcb
        if ~ifs
            wig = double(wig);
        elseif ifs==-1
            wig = simplify(wig);
        end
    end

    if ~isempty(NKs) % incorrect numerical results - swich to the symbolic computations
        NKs=sort(NKs);
        wig(NKs) = Wigner6j(j1(NKs),j2(NKs),j3(NKs),j4(NKs),j5(NKs),j6(NKs),0);
    end

    if ~isempty(kk)
        k0 = diff(kk);
        if all(k0~=1)
            wig(kk) = wig(kk+1);
        else
            k0 = find(k0>1);
            k00=1;
            for k=1:numel(k0)
                wig(kk(k00:k0(k))) = wig(kk(k0(k))+1);
                k00 = k0(k)+1;
            end
            wig(kk(k00:end)) = wig(kk(end)+1);
        end
     end

    if exist('ind0','var')
        wig(ind0)=wig;
    end

    if ifcb % recompute to the coupling matrix element
        if ifs>0
            wig = wig .* sqrt((2*J3+1).*(2*J6+1)) .* (-1).^(J1+J2+J4+J5);
        elseif ifs==0
            wig = double( wig .* sqrt(sym((2*J3+1).*(2*J6+1))) .* (-1).^(J1+J2+J4+J5) );
        elseif ifs==-1
            wig = simplify( wig .* sqrt(sym((2*J3+1).*(2*J6+1))) .* (-1).^(J1+J2+J4+J5) );
        else
            wig = wig .* sqrt(sym((2*J3+1).*(2*J6+1))) .* (-1).^(J1+J2+J4+J5);
        end
    end
    
    if NJ>1 && prod(Nj1S)==NJ
        wig=reshape(wig,Nj1S);
    end

catch mectc
%%
    beep;
    disp(mectc);
    disp('Wigner6j ERROR: The input arguments can be incorrect; see comments in the code,');
    if ifs<=0 || exist('NKs','var') && ~isempty(NKs)
        disp('or the ``Symbolic'' toolbox is not properly installed')
    end
    try
        if NJ>1 && prod(Nj1S)==NJ
            wig=NaN(Nj1S);
        else
            wig=NaN(size(j1));
        end
    catch
        wig=NaN;
    end
    beep;
end

return;
end
