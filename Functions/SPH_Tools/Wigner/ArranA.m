function [ind0,varargout] = ArranA (varargin)
% simultaneous sorting of several input column vectors
% so that varargout{1} becomes an acsending-order version of varargin{1},
% and varargout{2} is resorted so that elements of varargin{2} corresponding
% to equal elements of i1 becomes an ascending-order subvectors as well, etc.;
% ind0  is a permutation vector.
%
% USAGE: [ind0,i1,i2,...] = ArranA (j1,j2,...)
% with j1, j2, ... being numerical arrays (matrices) containing the vectors to be arranged column-wise.
%
% Programmer: V. B. Sovkov
% St. Petersburg State University
% Shanxi University
%%
try
    if ~nargin
        ind0=[];
        varargout={};
        return;
    end
    [varargout{nargin},ind0] = sort(varargin{nargin},1);
    for m=1:size(ind0,2)
        for k=nargin-1:-1:1
            varargout{k}(:,m) = varargin{k}(ind0(:,m),m);
        end
    end
    for n=nargin-1:-1:1
        [varargout{n},ind] = sort(varargout{n},1);
        for m=1:size(ind0,2)
            for k=nargin:-1:1
                if k~=n
                    varargout{k}(:,m) = varargout{k}(ind(:,m),m);
                end
            end
            ind0(:,m)=ind0(ind(:,m),m);
        end
    end
catch mectc
    varargout=[];
    disp(mectc);
end
return