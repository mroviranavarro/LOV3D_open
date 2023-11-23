function wig = Wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9,ifs,ifcb)
% Computes the Wigner 9j symbols
%    wig = {j1 j2 j3
%           j4 j5 j6
%           j7 j8 j9}
% or the recoupling matrix elements
%    wig = <(j1,j2)j3,(j4,j5)j6,j9|(j1,j5)j7,(j2,j4)j8,j9>
% using their expression in terms of the 6j-symbols.
%
%     USAGE
% wig = Wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9,ifs,ifcb);
% Wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
%
%     INPUT
%     mandatory arguments:
% j1, j2, j3, j4, j5, j6, j7, j8, j9
%      - the sets of the angular momenta; arrays of the same sizes;
%     optional arguments:
% ifs  - the switch of the computational methods;
%        the main use of the ifs parameter is its transferring to the 6j-symbol procedure---see comments in Wigner6j.m for details;
%        although the central part of all the computations is based on the Racah formula, some details differ:
%      = 0 - the symbolic (presumably accurate) computation with the double-precision output;
%      =-1 - the same symbolic computation as ifs=0 but with the symbolic output (accurate square root/rational-type equation, simplified);
%      =-2 - the same as ifs=-1 but without a final simplification of the symbolic expression, which can be time-consuming sometimesl
%      = 1, =2 (default, recommended), =3 - numeric double-precision algorithms (see REMARKS below);
%      all other input values of ifs are set to the closest from the above list.
% ifcb - if exists and is true, switches to computing the coupling matrix elements instead of the 9j-symbols;
%        default ifcb=false (9j-symbols are computed);
%     OUTPUT
% wig  - the resulting values of the 9j-symbols (ifcb=false) or the coupling matrix element (ifcb=true) in either numeric double-precision or symbolic form (see the input parameter ifs);
%        array of the same size as j1.
%
%     LINKS (DEFINITIONS)
% (1) https://en.wikipedia.org/wiki/9-j_symbol
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

    NJ = min([numel(j1) numel(j2) numel(j3) numel(j4) numel(j5) numel(j6) numel(j7) numel(j8) numel(j9)]);
    Nj1S=size(j1);
    j1 = double(reshape(j1(1:NJ),NJ,1));
    j2 = double(reshape(j2(1:NJ),NJ,1));
    j3 = double(reshape(j3(1:NJ),NJ,1));
    j4 = double(reshape(j4(1:NJ),NJ,1));
    j5 = double(reshape(j5(1:NJ),NJ,1));
    j6 = double(reshape(j6(1:NJ),NJ,1));
    j7 = double(reshape(j7(1:NJ),NJ,1));
    j8 = double(reshape(j8(1:NJ),NJ,1));
    j9 = double(reshape(j9(1:NJ),NJ,1));

    if nargin<10 || isempty(ifs) || ~isnumeric(ifs) && ~islogical(ifs) || ifs(1)>0 && ifs(1)<=1.5 % choose the algorithm
        ifs = 2; % default
        ifs0=ifs;
    elseif ifs(1)>0
        ifs = min(round(ifs(1)),3);
        ifs0=ifs;
    else
        ifs=max(round(ifs(1)),-2);
        ifs0=-2;
    end
    
    NK = if3jc(j1(1:NJ),j2(1:NJ),j3(1:NJ)) & if3jc(j4(1:NJ),j5(1:NJ),j6(1:NJ)) & if3jc(j7(1:NJ),j8(1:NJ),j9(1:NJ)) & if3jc(j1(1:NJ),j4(1:NJ),j7(1:NJ)) & if3jc(j2(1:NJ),j5(1:NJ),j8(1:NJ)) & if3jc(j3(1:NJ),j6(1:NJ),j9(1:NJ));
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

    ifcb = nargin>10 && ~isempty(ifcb) && ifcb(1); % if the Clebsch-Gordan coefficient instead of the 3j-symbols are needed?
    if ifcb
%         J1=j1;
%         J2=j2;
        J3=j3;
%         J4=j4;
%         J5=j5;
        J6=j6;
        J7=j7;
        J8=j8;
%         J9=j9;
    end

    s = zeros(NJ,1);
    for k=1:NJ
        ind1=[2,3,1];
        ind2=ind1;
        while ~issorted(ind1) || ~issorted(ind2)
            j01=[j1(k);j2(k);j3(k)];
            j02=[j4(k);j5(k);j6(k)];
            j03=[j7(k);j8(k);j9(k)];
            [ind1,j01,j02,j03] = ArranA(j01,j02,j03);
            if ~issorted(ind1)
                j1(k)=j01(1);
                j2(k)=j01(2);
                j3(k)=j01(3);
                j4(k)=j02(1);
                j5(k)=j02(2);
                j6(k)=j02(3);
                j7(k)=j03(1);
                j8(k)=j03(2);
                j9(k)=j03(3);
                s(k) = ceil(mod(s(k)+sper(ind1),2));
            end
            j01=[j1(k);j4(k);j7(k)];
            j02=[j2(k);j5(k);j8(k)];
            j03=[j3(k);j6(k);j9(k)];
            [ind2,j01,j02,j03] = ArranA(j01,j02,j03);
            if ~issorted(ind2)
                j1(k)=j01(1);
                j4(k)=j01(2);
                j7(k)=j01(3);
                j2(k)=j02(1);
                j5(k)=j02(2);
                j8(k)=j02(3);
                j3(k)=j03(1);
                j6(k)=j03(2);
                j9(k)=j03(3);
                s(k) = ceil(mod(s(k)+sper(ind2),2));
            end
        end
    end
    if NJ>1
        [ind0,j1,j2,j3,j4,j5,j6,j7,j8,j9] = ...
            ArranA ( reshape(j1(1:NJ),NJ,1),reshape(j2(1:NJ),NJ,1),reshape(j3(1:NJ),NJ,1),reshape(j4(1:NJ),NJ,1),reshape(j5(1:NJ),NJ,1),reshape(j6(1:NJ),NJ,1),reshape(j7(1:NJ),NJ,1),reshape(j8(1:NJ),NJ,1),reshape(j9(1:NJ),NJ,1));
        NK = NK(ind0);
        s = s(ind0);
    end
    sj = mod(j1+j2+j3+j4+j5+j6+j7+j8+j9,2);
    s = s & sj;

    kk = find ( sj & (j1==j2 & j4==j5 & j7==j8 | j2==j3 & j5==j6 & j8==j9 | j1==j4 & j2==j5 & j3==j6 | j4==j7 & j5==j8 & j6==j9) );
    if ~isempty(kk)
        NK(kk)=0;
    end

    kk = find ( j1(1:NJ-1)==j1(2:NJ) & j2(1:NJ-1)==j2(2:NJ) & j3(1:NJ-1)==j3(2:NJ) & j4(1:NJ-1)==j4(2:NJ) & j5(1:NJ-1)==j5(2:NJ) & j6(1:NJ-1)==j6(2:NJ) );
    if ~isempty(kk)
        NK(kk)=0;
    end

    NK = find(NK);

%%

    k0 = find(j1(NK)==0);
    if ~isempty(k0)
        k = NK(k0);
        k = k(j4(k)==j7(k) & j2(k)==j3(k));
        if ~isempty(k)
            if ifs>0
                wig(k) = Wigner6j(j2(k),j5(k),j8(k),j4(k),j9(k),j6(k),ifs0) ./...
                    sqrt((2*j4(k)+1).*(2*j2(k)+1)) .* (-1).^(j6(k)+j4(k)+j8(k)+j2(k)+s(k));
            else
                wig(k) = Wigner6j(j2(k),j5(k),j8(k),j4(k),j9(k),j6(k),ifs0) ./...
                    sqrt(sym((2*j4(k)+1).*(2*j2(k)+1))) .* (-1).^(j6(k)+j4(k)+j8(k)+j2(k)+s(k));
            end
        end
        NK(k0)=[];
    end
    if ~isempty(NK)
        t2 = min([j1(NK)+j9(NK),j2(NK)+j6(NK),j4(NK)+j8(NK)],[],2);
        t1 = max(abs([j1(NK)-j9(NK),j2(NK)-j6(NK),j4(NK)-j8(NK)]),[],2);
        for k=1:numel(NK)
            kN=NK(k);
            t = t1(k):t2(k);
            k12=numel(t);
            if k12
                w = Wigner6j(...
                    [ repmat(j1(kN),1,k12);       repmat(j2(kN),1,k12);       repmat(j3(kN),1,k12) ],...
                    [ repmat(j4(kN),1,k12);       repmat(j5(kN),1,k12);       repmat(j6(kN),1,k12) ],...
                    [ repmat(j7(kN),1,k12);       repmat(j8(kN),1,k12);       repmat(j9(kN),1,k12) ],...
                    [ repmat(j8(kN),1,k12);       repmat(j4(kN),1,k12);               t            ],...
                    [ repmat(j9(kN),1,k12);               t           ;       repmat(j1(kN),1,k12) ],...
                    [         t           ;       repmat(j6(kN),1,k12);       repmat(j2(kN),1,k12) ],...
                    ifs0);
                wig(kN) = prod(w,1) .* (-1).^(s(kN)+2*t) * (2*t'+1);
            end
        end
    end

%   post-processing

    if ~ifcb
        if ~ifs
            wig = double(wig);
        elseif ifs==-1
            wig = simplify(wig);
        end
    end

    if ~isempty(kk)
        k0 = diff(kk);
        if all(k0~=1)
            wig(kk) = (-1).^(s(kk)-s(kk+1)) .* wig(kk+1);
        else
            k0 = find(k0>1);
            k00=1;
            for k=1:numel(k0)
                wig(kk(k00:k0(k))) = (-1).^(s(kk(k00:k0(k)))-s(kk(k0(k))+1)) .* wig(kk(k0(k))+1);
                k00 = k0(k)+1;
            end
            wig(kk(k00:end)) = (-1).^(s(kk(k00:end))-s(kk(end)+1)) .* wig(kk(end)+1);
        end
     end

    if exist('ind0','var')
        wig(ind0)=wig;
    end

    if ifcb % recompute to the coupling matrix element
        if ifs>0
            wig = wig .* sqrt((2*J3+1).*(2*J6+1).*(2*J7+1).*(2*J8+1));
        elseif ifs==0
            wig = double( wig .* sqrt(sym((2*J3+1).*(2*J6+1).*(2*J7+1).*(2*J8+1)))  );
        elseif ifs==-1
            wig = simplify( wig .* sqrt(sym((2*J3+1).*(2*J6+1).*(2*J7+1).*(2*J8+1))) );
        else
            wig = wig .* sqrt(sym((2*J3+1).*(2*J6+1).*(2*J7+1).*(2*J8+1)));
        end
    end
    
    if ifs>0
        k=find(isinf(wig) | isnan(wig));
        if ~isempty(k)
            wig(k) = (-1).^s(k) .* Wigner9j(j1(k),j2(k),j3(k),j4(k),j5(k),j6(k),j7(k),j8(k),j9(k),0,ifcb);
        end
    end
    
    if NJ>1 && prod(Nj1S)==NJ
        wig=reshape(wig,Nj1S);
    end

catch mectc
%%
    beep;
    disp(mectc);
    disp('Wigner9j ERROR: The input arguments can be incorrect; see comments in the code,');
    if ifs<=0 % || exist('NKs','var') && ~isempty(NKs)
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
