% *************************************************************************
%
% ULTRA-SR | Check Localization Accuracy
%
% Inputs
%   truthPoints - True bubble locations [N-by-2] OR [N-by-3], but only in-
%                 plane (x-z plane) will be checked [m]
%   testPoints  - Localized bubble locations [N-by-2] OR [N-by-3] [m]
%   tolerance   - If the test is within this distance of truth, the
%                 localization is treated as correct.
%
% Outputs
%   correct  - 1-by-N boolean vector if peak is close enough to truth
%   distance - 1-by-N vector of Euclidean distances.
%   
% For further details, please see (and consider citing if used)
%   [1] "MR for ULTRA-SR: Improved Localization with Morphological Image 
%       Processing". IEEE IUS Proc. (2022) 
%       DOI: 10.1109/IUS54386.2022.9957276
%
%         Scott Schoen Jr | MGH-CURT | sschoenjr@mgh.harvard.edu
%
% *************************************************************************

function [correct, distances] = checkProximity(truthPoints, testPoints, tolerance)

% First, get distances
[closestInds, distances] = dsearchn(truthPoints, testPoints);

% Keep only the closest point if two test points are matched to the same
% truth
for indCount = 1 : length(closestInds)
    
    % If we've already removed this point, go on
    if isnan(closestInds(indCount))
        continue;
    end
    
    % Find all found peaks matched to this truth
    repeatInds = find( closestInds == closestInds(indCount) );
    
    % Keep only the closest one
    [~, minInd] = min( closestInds(repeatInds) );
    repeatInds( minInd ) = [];
    closestInds(repeatInds) = NaN;
    
end

% Return true/false vector
correct = 0.*distances;
for cCount = 1 : length( closestInds )
    
    % If the closest truth point was already used, automatically false.
    if isnan( closestInds(cCount) )
        continue;
    end
    
    % Get index of closest point in truth set
    loopInd = closestInds(cCount);
    
    % Get distance from test point to closest point
    r0 = [ truthPoints( loopInd, 1 ), truthPoints( loopInd, 2 ) ];
    r1 = [ testPoints( cCount, 1 ), testPoints( cCount, 2 ) ];
    distToPoint = sqrt( sum( (r1 - r0).^(2) ) );
    
    % Check if it's close enough
    pointCloseEnough = distToPoint < tolerance;
    if pointCloseEnough
        correct( cCount ) = 1;
    end
    
end

end

