function tri = tc(j1,j2,j3,ifs)
% sqrt of the triangle coefficient or its logarithm (when ifs=3)
if nargin<4 || isempty(ifs) || ~isnumeric(ifs) || ~isreal(ifs) || ifs>0
    if ifs==3
        tri = zeros(size(j1));
        for k=1:numel(tri)
            tri(k) = gammaln(j1(k)+j2(k)-j3(k)+1)/2 + gammaln(j1(k)-j2(k)+j3(k)+1)/2 - sum( log( (-j1(k)+j2(k)+j3(k)+1):(j1(k)+j2(k)+j3(k)+1) )/2 );
        end
    else
        if ifs==1
            tri = sqrt( factorial(j1+j2+j3+1) ./ factorial(j1+j2-j3) );
        elseif ifs==2
            for k=numel(j1):-1:1
                tri(k,1) = prod( sqrt( (j1(k)+j2(k)-j3(k)+1):(j1(k)+j2(k)+j3(k)+1) ) );
            end
        end
        k=find(~isinf(tri) & ~isnan(tri));
        if ~isempty(k)
            tri(k) = sqrt( factorial(-j1(k)+j2(k)+j3(k)) ) ./ tri(k) .* sqrt( factorial(j1(k)-j2(k)+j3(k)) );
        end
    end
else
    tri = sqrt( factorial(sym(j1+j2-j3))./ factorial(sym(j1+j2+j3+1)) .* factorial(sym(j1-j2+j3)) .* factorial(sym(-j1+j2+j3)) );
end
return;
end
