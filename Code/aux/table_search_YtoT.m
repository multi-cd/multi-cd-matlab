function Tsol = table_search_YtoT(Ytarget,Ygrid,Tgrid,varargin)
% table_search_YtoT: Performs a simple table search with interpolation.
% Assume that there is a continuous (and monotonic) relationship
% between two variables T and Y. 
% From a pre-computed table for (T,Y) values that samples the relationship,
% we find Tsol that corresponds to the requested Ytarget, as (Tsol,Ytarget).
% 
% INPUT:
%   - Ytarget: the target value of Y (single number)
%   - Ygrid: pre-computed values of Y (to be paired with Tgrid)
%   - Tgrid: pre-computed values of T (to be paired with Ygrid)
%   - myTmax: [OPTIONAL] fixed value for T if Ytarget exceeds Ygrid range
% 
% OUTPUT:
%   - Tsol: the "solved" value of T, such that (Tsol,Ytarget) is a pair
% ------------------------------------------------------------------------

% 2018 Ji Hyun Bak


%% detect data size

dataSizeCutoff = 1000; % some big number (1000 is a small cutoff)
if(numel(Ytarget)>dataSizeCutoff)
    method = 'largedata';
else
    method = 'smalldata';
end

%% table search and interpolation

% locate target bin in the grid
[ilower,iupper] = locateIdxFromTable(Ytarget,Ygrid,method);

% locate target point fraction within the Y bin (interpolation)
Ylower = Ygrid(ilower);
Yupper = Ygrid(iupper);
Tfrac = (Ytarget(:)-Ylower)./(Yupper-Ylower); 

% project back to T space
Tlower = Tgrid(ilower);
Tupper = Tgrid(iupper);
Tsol = Tlower + (Tupper-Tlower).*Tfrac;

% finally, fix where Y was beyond the max of tabular value 
if(nargin>3)
    myTmax = varargin{1}; % optional input to specify Tmax
else
    % if not specified, use default value
    % myTmax = Inf;
    myTmax = max(Tgrid); 
end
Tsol(Ytarget>=max(Ygrid)) = myTmax;   % bug fix 2/21/2019 (replaced > with >=)

end

% ------------------------------------------------------------------------

function [ilower,iupper] = locateIdxFromTable(Ytarget,Ygrid,method)
% locate target bin in the grid

ngrid = numel(Ygrid);

switch lower(method)
    case 'smalldata'
        
        compmat = bsxfun(@ge,Ytarget(:),Ygrid(:)'); % compare targets to Ygrid table
        icomp = sum(compmat,2); % can range from 0 to ngrid
    case 'largedata'
        
        ntarg = numel(Ytarget);
        icomp = NaN(ntarg,1);
        for ii = 1:ntarg
            icomp(ii) = sum(Ytarget(ii)>=Ygrid); % can range from 0 to ngrid
        end
    otherwise
        error('table search: unknown method');
end

ilower = max(icomp,1); % index of lower bin
iupper = min(icomp+1,ngrid); % index of upper bin

end
