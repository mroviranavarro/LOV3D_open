function if3j = if3jc(J1,J2,J3)
% checks necessary conditions on the angular momenta J1, J2, J3, which must be:
% (1) non-negative;
% (2) either integer or half-integer;
% (3) their sum must be integer;
% (4) triangular rule.
% returns true if all the conditions are fulfilled; otherwise returns false.

if3j = J1(:)>=0 & J2(:)>=0 & J3(:)>=0 & mod(J1(:)-J2(:)-J3(:),1)==0 & mod(2*J1(:),1)==0 & mod(2*J2(:),1)==0 & mod(2*J3(:),1)==0 & J3(:) <= J1(:)+J2(:) & J3(:) >= abs(J1(:)-J2(:));

return;
end


