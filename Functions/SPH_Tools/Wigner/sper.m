function s=sper(ind)
% permutation parity
if size(ind,1)==1
    ind=ind';
end
s = zeros(1,size(ind,2));

for k=1:size(ind,1)-1
    s = s + sum(ind(k,:) > ind(k+1:end,:),1);
end
s = mod(s,2);
return;
end
