% *************************************************************************
%
% ULTRA-SR | Create Accumulation Image
%
%   Function reads in SR peaks and assigns a small Gaussian distirbution to
%   each to create a super-resolution accumulation map.
%
% Inputs
%   xPeaks   - Detected x position of peaks [m]
%   zPeaks   - Detected z position of peaks [m]
%   sz       - Standard deviation(s) of the Gaussian PSFs to use [m, (m)]
%   x        - x grid to define final image over [m]
%   z        - z grid to define final image over [m]
%   ditherPoints
%            - If positive, a small random displacement is added to each to
%              reduce gridding artifact [0, 0.5]
%   ampPeaks - Whether or not to weight peaks by amplitude. Recommend
%              false.
%
% Outputs
%   img - SR accumulation map.
%   
% For further details, please see (and consider citing if used)
%   [1] "MR for ULTRA-SR: Improved Localization with Morphological Image 
%       Processing". IEEE IUS Proc. (2022) 
%       DOI: 10.1109/IUS54386.2022.9957276
%
%         Scott Schoen Jr | MGH-CURT | sschoenjr@mgh.harvard.edu
%
% *************************************************************************

function [img] = ...
    makeSrImg( xPeaks, zPeaks, sz, x, z, ditherPoints, ampPeaks )

% Check if we want to slightly dither points
if nargin < 6 || ~ditherPoints
    dispVal = 0;
else
    if isa( ditherPoints, 'double' )
        dispVal = min( 1, max( ditherPoints, 0 ) );
    else
        dispVal = 0.2;
    end
end


% Check if we want to apply weights
if nargin < 7
    peakVals = 0.*xPeaks + 1;
else
    if isa( ditherPoints, 'double' )
        peakVals = ampPeaks;
    else
        peakVals = 0.*xPeaks + 1;
    end
end

if numel(sz) == 2
    szx = sz(1);
    szz = sz(2);
else
    szx = sz;
    szz = sz;
end
    
% Get pixel size
dx = x(1,2) - x(1,1);
dz = z(2,1) - z(1,1);
dr = sqrt( dx.^(2) + dz.^(2) );

% Define PSF
gaussianDist = @(xc, zc) ...
    exp( -( (x - xc).^(2)./szx.^(2) + (z - zc).^(2)./szz.^(2) ) );

% Initialize
img = 0.*x;

% Plot each
numPeaks = length( xPeaks );
parfor peakCount = 1 : numPeaks    
    
    % Add small displacement if required
    r = dispVal.*rand(1).*dr;
    theta = 2.*pi.*rand(1);
    
    xC = xPeaks(peakCount) + r.*cos(theta);
    zC = zPeaks(peakCount) + r.*sin(theta);
    
    % Add contribution
    img = img + peakVals(peakCount).*gaussianDist( xC, zC );
    
end

end

