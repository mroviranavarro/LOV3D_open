function wig = Wigner3j(j1,j2,j3,m1,m2,m3,ifs,ifcb)
% Computes the Wigner 3j symbols
%    wig =(j1 j2 j3
%          m1 m2 m3)
% or the Clebsch-Gordan coefficient
%    wig = <j1,m1,j2,m2|j3,m3>
% using the Racah formula.
%
%     USAGE
% wig = Wigner3j(j1,j2,j3,m1,m2,m3,ifs,ifcb);
% Wigner3j(j1,j2,j3,m1,m2,m3)
%
%     INPUT
%     mandatory arguments:
% j1, j2, j3, m1, m2, m3
%      - the sets of the angular momenta and their projections; arrays of the same sizes;
%     optional arguments, the default values are adopted when they are absent or empty in the input:
% ifs  - the switch of the computational methods;
%        although the central part of all the computations is based on the Racah formula, some details differ:
%      = 0 - the symbolic (presumably accurate) computation with the double-precision output;
%      =-1 - the same symbolic computation as ifs=0 but with the symbolic output (accurate square root/rational-type equation, simplified);
%      =-2 - the same as ifs=-1 but without a final simplification of the symbolic expression, which can be time-consuming sometimesl
%      = 1 (default, recommended), =2, =3 - numeric double-precision algorithms (see REMARKS below);
%      all other input values of ifs are set to the closest from the above list.
% ifcb - if exists and is true, switches to computing the Clebsch-Gordan coefficients instead of the 3j-symbols;
%        default ifcb=false (3j-symbols are computed);
%     OUTPUT
% wig  - the resulting values of the 3j-symbols (ifcb=false) or the Clebsch-Gordan coefficients (ifcb=true) in either numeric double-precision or symbolic form (see the input parameter ifs);
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
% (g) the algorithm ifs=1 ensures a reasonably good accuracy; the worst registered case was
%     Wigner3j(56,59,51,-2,14,-12,1)=-0.00019235064 (exact Wigner3j(56,59,51,-2,14,-12,0)=-0.00019235487)
%     with the relative inaccuracy of 2.2e-5
%     (occuring before switching to the symbolic regime);
% (h) the algorithm ifs=2 (comparing to ifs=1) proceeds in a wider range of quantum numbers numerically
%     before switching to the symbolic regime (i.e., can work faster with big quantum numbers)
%     but with somewhat lower precision; in our experimental probing we registered the cases with the worst
%     relative inaccuracy of 0.13 for Wigner3j(80.5,85,68.5,-18.5,-15,33.5,2)=-2.781e-06 (exact Wigner3j(80.5,85,68.5,-18.5,-15,33.5,0)=-2.4637e-06)
%     and abolute inaccuracy of 0.0001 for Wigner3j(82.5,67,90.5,4.5,1,-5.5,2)=0.0045102 (exact 0.0046133=0.0046133);
%     in the range of the numeric calculations of the algorithm ifs=1, the algorithm ifs=2 provides the same characteristic accuracy as ifs=1;
% (i) the algorithm ifs=3 (comparing to ifs=2) proceeds in a wider range of quantum numbers numerically
%     before switching to the symbolic regime (i.e., can work faster with big quantum numbers)
%     but for some combinations of big quantum numbers, completely wrong results are observed:
%     e.g., Wigner3j(465.5,488.5,498,13.5,-90.5,77,3); Wigner3j(95,96,99,-17,1,16,3);
%     thereby, we do not recommend using the algorithm ifs=3 at all, and only keep it here in view of the completeness.
%
%     LINKS (DEFINITIONS)
% (1) https://en.wikipedia.org/wiki/3-j_symbol
% (2) https://en.wikipedia.org/wiki/Clebsch萌ordan_coefficients
% (3) R. N. Zare, Angular Momentum: Understanding Spatial Aspects in Chemistry and
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

%%
%   general check, preprocessing
    NJ = min([numel(j1) numel(j2) numel(j3) numel(m1) numel(m2) numel(m3)]);
    Nj1S=size(j1);
    j1 = double(reshape(j1(1:NJ),NJ,1));
    j2 = double(reshape(j2(1:NJ),NJ,1));
    j3 = double(reshape(j3(1:NJ),NJ,1));
    m1 = double(reshape(m1(1:NJ),NJ,1));
    m2 = double(reshape(m2(1:NJ),NJ,1));
    m3 = double(reshape(m3(1:NJ),NJ,1));
    
    if nargin<7 || isempty(ifs) || ~isnumeric(ifs) && ~islogical(ifs) || ifs(1)>0 && ifs(1)<=1.5 % choose the algorithm
        ifs = 1; % default
    elseif ifs(1)>0
        ifs = min(round(ifs(1)),3);
    else
        ifs=max(round(ifs(1)),-2);
    end

    ifcb = nargin>7 && ~isempty(ifcb) && ifcb(1); % if the Clebsch-Gordan coefficient instead of the 3j-symbols are needed?
    if ifcb
        m3 = -m3;
        M3=m3;
        J1=j1;
        J2=j2;
        J3=j3;
    end

    NK =  if3jc(j1,j2,j3) & (m1+m2+m3==0) & ( j1 - m1 == floor ( j1 - m1 ) ) & ( j2 - m2 == floor ( j2 - m2 ) ) & ( j3 - m3 == floor ( j3 - m3 ) ) & (abs(m1) <= j1) & (abs(m2) <= j2) & (abs(m3) <= j3);
    if all(~NK) % all the results are zero based on the selection rules
        if ifs>=0
            wig=zeros(NJ,1);
        else
            wig=sym(zeros(NJ,1));
        end
        return;
    end

    if ifs>0
        wig=zeros(NJ,1);
    else
        wig=sym(zeros(NJ,1));
    end

    sp = ones(NJ,1); % to keep the signs produced by the permutatiions

    for k=1:NJ % permute the elements of every set of quantum numbers into the standard (ascending) order
        jj = [j1(k) ; j2(k) ; j3(k)];
        mm = [m1(k) ; m2(k) ; m3(k)];
        [ind,a,b] = ArranA(jj,abs(mm));
        j1(k) = a(1);
        j2(k) = a(2);
        j3(k) = a(3);
        m1(k) = b(1);
        if mm(ind(1))<0 || (  ~mm(ind(1)) && ( mm(ind(2))<0 || (~mm(ind(2)) && mm(ind(3))<0) )  )
            m2(k) = -mm(ind(2));
            m3(k) = -mm(ind(3));
            if mod(j1(k)+j2(k)+j3(k),2) && ~sper(ind)
                sp(k) = -1;
            end
        else
            m2(k) = mm(ind(2));
            m3(k) = mm(ind(3));
            if mod(j1(k)+j2(k)+j3(k),2) && sper(ind)
                sp(k) = -1;
            end
        end
    end
    
    if NJ>1 % arrange the sets into the standard (ascending) order
        [ind0,j1,j2,j3,m1,m2,m3] = ...
            ArranA ( j1,j2,j3,m1,m2,m3);
        NK = NK(ind0);
        sp = sp(ind0);
    end

    kk = find ( mod(j1+j2+j3,2) & (m1==0 & m2==0 & m3==0 | j1==j2 & m1==m2 | j2==j3 & m2==m3) ); % find some obviously zero (due to symmetry properties) coefficients
    if ~isempty(kk)
        NK(kk)=0; % exclude the obviously zero coefficients
    end

    kk = find ( j1(1:NJ-1)==j1(2:NJ) & j2(1:NJ-1)==j2(2:NJ) & j3(1:NJ-1)==j3(2:NJ) & m1(1:NJ-1)==m1(2:NJ) & m2(1:NJ-1)==m2(2:NJ) & m3(1:NJ-1)==m3(2:NJ) ); % find duplicates
    if ~isempty(kk)
        NK(kk)=0; % exclude duplicates in order to avoid repetitive computations
    end

    NKs=[]; % for pointing to erronious numeric results in order to switch them to the symbolic computations
%%
%   special cases - faster and more accurate computation than the general case below

    if any(NK & j1<=2) % the range of the special cases currently included

%   case j = 0
        kp = find(NK & j1==0);
        if ~isempty(kp)
            if ifs>0
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp)) ./ sqrt(2*j2(kp)+1);
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp)) ./ sqrt(sym(2*j2(kp)+1));
            end
            NK(kp)=0;
        end
%   case   j = 1/2
        kp = find(NK & j1==1/2);
        if ~isempty(kp)
            if ifs>0
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)-1/2) .* sqrt( (j2(kp)-m3(kp)+0.5) ./ (2*j2(kp)+2) ./ (2*j2(kp)+1) );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)-1/2) .* sqrt( sym((j2(kp)-m3(kp)+1/2) ./ (2*j2(kp)+2) ./ (2*j2(kp)+1)) );
            end
            NK(kp)=0;
        end
%   case   j = 1
        kp = find(NK & j1==1 & j2==j3 & m1==0);
        if ~isempty(kp)
            if ifs>0
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp)) .* (2*m2(kp)) ./ sqrt( (2*j2(kp)+2) .* (2*j2(kp)+1)  .* (2*j2(kp)) );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp)) .* sym(2*m2(kp)) ./ sqrt( sym( (2*j2(kp)+2) .* (2*j2(kp)+1)  .* (2*j2(kp)) ) );
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & m1==0); %  & j2~=j3
        if ~isempty(kp)
            if ifs>0
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp)-1) .* sqrt( 2*(j2(kp)+m2(kp)+1) .* (j2(kp)-m2(kp)+1) ./ (2*j2(kp)+3) ./ (2*j2(kp)+2)  ./ (2*j2(kp)+1) );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp)-1) .* sqrt( sym ( 2*(j2(kp)+m2(kp)+1) .* (j2(kp)-m2(kp)+1) ./ (2*j2(kp)+3) ./ (2*j2(kp)+2)  ./ (2*j2(kp)+1) ) );
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & j2==j3 & m1==1);
        if ~isempty(kp)
            if ifs>0
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp)) .*  sqrt( 2*(j2(kp)-m2(kp)) .* (j2(kp)+m2(kp)+1) ./ (2*j2(kp)+2) ./ (2*j2(kp)+1)  ./ (2*j2(kp)) );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp)) .*  sqrt( sym( 2*(j2(kp)-m2(kp)) .* (j2(kp)+m2(kp)+1) ./ (2*j2(kp)+2) ./ (2*j2(kp)+1)  ./ (2*j2(kp)) ) );
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==1 & m1==1); %  & j2~=j3
        if ~isempty(kp)
            if ifs>0
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)-1) .* sqrt( (j2(kp)-m3(kp)) .* (j2(kp)-m3(kp)+1) ./ (2*j2(kp)+3) ./ (2*j2(kp)+2)  ./ (2*j2(kp)+1) );
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)-1) .* sqrt( sym ( (j2(kp)-m3(kp)) .* (j2(kp)-m3(kp)+1) ./ (2*j2(kp)+3) ./ (2*j2(kp)+2)  ./ (2*j2(kp)+1) ) );
            end
            NK(kp)=0;
        end
%   case   j = 3/2
        kp = find(NK & j1==3/2 & j3-j2==1/2 & m1==1/2);
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+3;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)+1/2)...
                    .* (j2(kp)+3*m3(kp)+3/2) .* sqrt(j2(kp)-m3(kp)+1/2)...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+3);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)+1/2)...
                    .* sym(j2(kp)+3*m3(kp)+3/2) .* sqrt(sym(j2(kp)-m3(kp)+1/2))...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & m1==1/2); %  & j3-j2==3/2
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+4;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)+1/2)...
                    .* sqrt( 3 * (j2(kp)-m3(kp)+1/2) .* (j2(kp)-m3(kp)+3/2) .* (j2(kp)+m3(kp)+3/2) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+4);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)+1/2)...
                    .* sqrt( sym( 3 * (j2(kp)-m3(kp)+1/2) .* (j2(kp)-m3(kp)+3/2) .* (j2(kp)+m3(kp)+3/2) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & j3-j2==1/2 & m1==3/2);
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+3;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)+1/2)...
                    .* sqrt( 3 * (j2(kp)-m3(kp)-1/2) .* (j2(kp)-m3(kp)+1/2) .* (j2(kp)+m3(kp)+3/2) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+3);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)+1/2)...
                    .* sqrt( sym( 3 * (j2(kp)-m3(kp)-1/2) .* (j2(kp)-m3(kp)+1/2) .* (j2(kp)+m3(kp)+3/2) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==3/2 & m1==3/2); %  & j3-j2==3/2
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+4;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)+1/2)...
                    .* sqrt( (j2(kp)-m3(kp)-1/2) .* (j2(kp)-m3(kp)+1/2) .* (j2(kp)-m3(kp)+3/2) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+4);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp)+1/2)...
                    .* sqrt( sym( (j2(kp)-m3(kp)-1/2) .* (j2(kp)-m3(kp)+1/2) .* (j2(kp)-m3(kp)+3/2) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) )...
                    ;
            end
            NK(kp)=0;
        end
%   case   j = 2, m = 0
        kp = find(NK & j1==2 & j3==j2 & m1==0);
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+3;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp))...
                    .* 2*(3*m2(kp).^2-j2(kp).*(j2(kp)+1))...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+3);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp))...
                    .* sym( 2*(3*m2(kp).^2-j2(kp).*(j2(kp)+1)) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & m1==0);
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+4;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp))...
                    .* 2*m3(kp) .* sqrt( 6 * (j2(kp)+m2(kp)+1) .* (j2(kp)-m2(kp)+1) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+4);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp))...
                    .* sym(2*m3(kp)) .* sqrt( sym( 6 * (j2(kp)+m2(kp)+1) .* (j2(kp)-m2(kp)+1) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & m1==0); %  & j3-j2==2
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+5;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sqrt( 6 * (j2(kp)+m3(kp)+2) .* (j2(kp)-m3(kp)+2) .* (j2(kp)+m3(kp)+1) .* (j2(kp)-m3(kp)+1) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+5);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sqrt( sym( 6 * (j2(kp)+m3(kp)+2) .* (j2(kp)-m3(kp)+2) .* (j2(kp)+m3(kp)+1) .* (j2(kp)-m3(kp)+1) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end
%   case   j = 2, m = 1
        kp = find(NK & j1==2 & j3==j2 & m1==1);
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+3;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp))...
                    .* (2*m2(kp)+1) .* sqrt( 6 * (j2(kp)+m2(kp)+1) .* (j2(kp)-m2(kp)) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+3);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp))...
                    .* sym(2*m2(kp)+1) .* sqrt( sym( 6 * (j2(kp)+m2(kp)+1) .* (j2(kp)-m2(kp)) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & m1==1);
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+4;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* 2*(j2(kp)+2*m3(kp)+2) .* sqrt( (j2(kp)-m3(kp)+1) .* (j2(kp)-m3(kp)) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+4);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sym(2*(j2(kp)+2*m3(kp)+2)) .* sqrt( sym( (j2(kp)-m3(kp)+1) .* (j2(kp)-m3(kp)) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & m1==1); %  & j3-j2==2
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+5;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sqrt( 4 * (j2(kp)+m3(kp)+2) .* (j2(kp)-m3(kp)+2) .* (j2(kp)-m3(kp)+1) .* (j2(kp)-m3(kp)) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+5);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sqrt( sym( 4 * (j2(kp)+m3(kp)+2) .* (j2(kp)-m3(kp)+2) .* (j2(kp)-m3(kp)+1) .* (j2(kp)-m3(kp)) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end
%   case   j = 2, m = 2
        kp = find(NK & j1==2 & j3==j2 & m1==2);
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+3;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp))...
                    .* sqrt( 6 * (j2(kp)-m2(kp)-1) .* (j2(kp)-m2(kp)) .* (j2(kp)+m2(kp)+1) .* (j2(kp)+m2(kp)+2) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+3);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)-m2(kp))...
                    .* sqrt( sym( 6 * (j2(kp)-m2(kp)-1) .* (j2(kp)-m2(kp)) .* (j2(kp)+m2(kp)+1) .* (j2(kp)+m2(kp)+2) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & j3-j2==1 & m1==2);
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+4;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sqrt( 4 * (j2(kp)-m3(kp)-1) .* (j2(kp)-m3(kp)) .* (j2(kp)-m3(kp)+1) .* (j2(kp)+m3(kp)+2) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+4);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sqrt( sym( 4 * (j2(kp)-m3(kp)-1) .* (j2(kp)-m3(kp)) .* (j2(kp)-m3(kp)+1) .* (j2(kp)+m3(kp)+2) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end
        kp = find(NK & j1==2 & m1==2); %  & j3-j2==2
        if ~isempty(kp)
            if ifs>0
                A = 2*j2(kp)+5;
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sqrt( (j2(kp)-m3(kp)-1) .* (j2(kp)-m3(kp)) .* (j2(kp)-m3(kp)+1) .* (j2(kp)-m3(kp)+2) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
                k0=find(isinf(wig(kp)) | isnan(wig(kp)));
                if ~isempty(k0)
                    NKs = [NKs;NK(kp(k0))];
                end
            else
                A = sym(2*j2(kp)+5);
                wig(kp) = sp(kp) .* (-1).^(j2(kp)+m3(kp))...
                    .* sqrt( sym( (j2(kp)-m3(kp)-1) .* (j2(kp)-m3(kp)) .* (j2(kp)-m3(kp)+1) .* (j2(kp)-m3(kp)+2) ) )...
                    ./ sqrt( A .* (A-1) .* (A-2) .* (A-3) .* (A-4) )...
                    ;
            end
            NK(kp)=0;
        end

    end
%   If other exact equation are known for some other partial cases, they can be included in the same manner

%%
%   main computation --- general case

    NK = find(NK);

    if ~isempty(NK)

        if ifs>0 % numeric computations
            C=zeros(length(NK),1);
            if ifs==3
                for k=length(NK):-1:1
                    C(k) = sum(log(max(2,(-j1(NK(k))+j2(NK(k))+j3(NK(k))+1)):(j1(NK(k))+j2(NK(k))+j3(NK(k))+1)));
                    if isinf(C(k)) || isnan(C(k))
                        NKs = sort([NKs;NK(k)]);
                        NK(k) = [];
                        C(k)=[];
                    else
                        C(k)= (...
                              2*sum(  log( 2:( j1(NK(k))-abs(m1(NK(k))) ) )  )...
                            + 2*sum(  log( 2:( j2(NK(k))-abs(m2(NK(k))) ) )  )...
                            + 2*sum(  log( 2:( j3(NK(k))-abs(m3(NK(k))) ) )  )...
                            + sum(  log( (j1(NK(k))-abs(m1(NK(k)))+1):(j1(NK(k))+abs(m1(NK(k)))) )  )...
                            + sum(  log( (j2(NK(k))-abs(m2(NK(k)))+1):(j2(NK(k))+abs(m2(NK(k)))) )  )...
                            + sum(  log( (j3(NK(k))-abs(m3(NK(k)))+1):(j3(NK(k))+abs(m3(NK(k)))) )  )...
                            - C(k) + (...
                              sum(  log( 2:(j1(NK(k))+j2(NK(k))-j3(NK(k))) )  ) + sum(  log( 2:(j1(NK(k))-j2(NK(k))+j3(NK(k))) )  )...
                              ) )/2;
                        if isinf(C(k)) || isnan(C(k))
                            NKs = sort([NKs;NK(k)]);
                            NK(k) = [];
                            C(k)=[];
                        else
                        end
                    end
                end
                k0 = find(isinf(C) | isnan(C));
                if ~isempty(k0)
                    NKs = sort([NKs;NK(k0)]);
                    NK(k0) = [];
                    C(k0)=[];
                end
            else
                if ifs==2
                    for k=1:length(NK)
                        C(k) = prod((-j1(NK(k))+j2(NK(k))+j3(NK(k))+1):(j1(NK(k))+j2(NK(k))+j3(NK(k))+1));
                    end
                    k0 = find(isinf(C) | isnan(C));
                    if ~isempty(k0)
                        NKs = sort([NKs;NK(k0)]);
                        NK(k0) = [];
                        C(k0)=[];
                    end
                else % ifs==1 --- the default recommended case
                    C = factorial(j1(NK)+j2(NK)+j3(NK)+1);
                    k0 = find(isinf(C) | isnan(C));
                    if ~isempty(k0)
                        NKs = sort([NKs;NK(k0)]);
                        NK(k0) = [];
                        C(k0)=[];
                    end
                    if ~isempty(NK)
                        C = C ./ factorial(-j1(NK)+j2(NK)+j3(NK));
                    end
                end
                if ~isempty(NK)
                    C = (-1).^(j1(NK)-j2(NK)-m3(NK)) .* sp(NK) ...
                            .* sqrt( factorial(j1(NK)-j2(NK)+j3(NK)) ./ C )...
                            .* sqrt( factorial(j1(NK)+m1(NK)) ) .* sqrt( factorial(j1(NK)-m1(NK)) )...
                            .* sqrt( factorial(j2(NK)+m2(NK)) ) .* sqrt( factorial(j2(NK)-m2(NK)) )...
                            .* sqrt( factorial(j3(NK)+m3(NK)) ) .* sqrt( factorial(j3(NK)-m3(NK)) )...
                            .* sqrt( factorial(j1(NK)+j2(NK)-j3(NK)) )...
                            ;
                    k0=find(isinf(C) | isnan(C));
                    if ~isempty(k0)
                        C(k0)=[];
                        NKs = sort([NKs;NK(k0)]);
                        NK(k0)=[];
                    end
                end
            end
        else % symbolic computations
            C = (-1).^(j1(NK)-j2(NK)-m3(NK)) .* sp(NK) .* sqrt (  ...
                       factorial(sym(-j1(NK)+j2(NK)+j3(NK)))...
                    ./ factorial(sym(j1(NK)+j2(NK)+j3(NK)+1))...
                    .* factorial(sym(j1(NK)-j2(NK)+j3(NK)))...
                    .* factorial(sym(j1(NK)+j2(NK)-j3(NK)))...
                    .* factorial(sym(j3(NK)+m3(NK)))...
                    .* factorial(sym(j3(NK)-m3(NK)))...
                    .* factorial(sym(j2(NK)+m2(NK)))...
                    .* factorial(sym(j2(NK)-m2(NK)))...
                    .* factorial(sym(j1(NK)+m1(NK)))...
                    .* factorial(sym(j1(NK)-m1(NK)))...
                    );
        end

        if ~isempty(NK) % could be emptied due to Infs in the double-precision mode
            t1 = j2(NK) - m1(NK) - j3(NK);
            t2 = j1(NK) + m2(NK) - j3(NK);
            t3 = j1(NK) + j2(NK) - j3(NK);
            t4 = j1(NK) - m1(NK);
            t5 = j2(NK) + m2(NK);

            tmin = max( [zeros(size(t1)), t1, t2], [], 2 );
            tmax = min( [t3, t4, t5], [], 2 );

            for k=numel(NK):-1:1 % t below can be of different lengths for the different k! Cannot apply the element-wise operations in this sense!

                kN = NK(k);

                t = tmin(k):tmax(k);

                if ifs>0
                    if ifs==3
                        wig(kN) = sp(kN) * (-1).^(j1(kN)-j2(kN)-m3(kN)+t) * exp(C(k) - (...
                              gammaln(t-t1(k)+1) + gammaln(t-t2(k)+1) + gammaln(t+1)...
                            + gammaln(t3(k)-t+1) + gammaln(t4(k)-t+1) + gammaln(t5(k)-t+1)...
                            ))';
                    else
                        wig(kN) = (-1).^t * ( C(k) ...
                            ./ factorial(t)...
                            ./ factorial(t-t1(k))...
                            ./ factorial(t-t2(k))...
                            ./ factorial(t3(k)-t)...
                            ./ factorial(t4(k)-t)...
                            ./ factorial(t5(k)-t)...
                            )';
                    end
                else
                    wig(kN) = (-1).^t * ( C(k) ...
                        ./ factorial(sym(t))...
                        ./ factorial(sym(t-t1(k)))...
                        ./ factorial(sym(t-t2(k)))...
                        ./ factorial(sym(t3(k)-t))...
                        ./ factorial(sym(t4(k)-t))...
                        ./ factorial(sym(t5(k)-t))...
                        )';
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
        wig(NKs) = sp(NKs) .* Wigner3j(j1(NKs),j2(NKs),j3(NKs),m1(NKs),m2(NKs),m3(NKs),0);
    end

    if ~isempty(kk) % duplicates
        k0 = diff(kk);
        if all(k0~=1)
            wig(kk) = sp(kk+1).*sp(kk).*wig(kk+1);
        else
            k0 = find(k0>1);
            k00=1;
            for k=1:numel(k0)
                wig(kk(k00:k0(k))) = sp(kk(k0(k))+1).*sp(kk(k00:k0(k))).*wig(kk(k0(k))+1);
                k00 = k0(k)+1;
            end
            wig(kk(k00:end)) = sp(kk(end)+1).*sp(kk(k00:end)).*wig(kk(end)+1);
        end
    end

    if exist('ind0','var') % rearrange to the initial order
        wig(ind0)=wig;
    end

    if ifcb % recompute to the Clebsch-Gordan coefficients
        if ifs>0
            wig = wig .* sqrt(2*J3+1) .* (-1).^(J1-J2-M3);
        elseif ifs==0
            wig = double( wig .* sqrt(sym(2*J3+1)) .* (-1).^(J1-J2-M3) );
        elseif ifs==-1
            wig = simplify(wig .* sqrt(sym(2*J3+1)) .* (-1).^(J1-J2-M3));
        else
            wig = wig .* sqrt(sym(2*J3+1)) .* (-1).^(J1-J2-M3);
        end
    end
    
    if NJ>1 && prod(Nj1S)==NJ
        wig=reshape(wig,Nj1S);
    end

catch mectc
%%
    beep;
    disp(mectc);
    disp('Wigner3j ERROR: The input arguments can be incorrect; see comments in the code,');
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

