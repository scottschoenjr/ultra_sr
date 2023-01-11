% *************************************************************************
%
% ULTRA-SR | Bubble Localization
%
% Inputs
%   x   - Meshgrid of x-positions in the image [m]
%   z   - Meshgrid of z-positions in the image [m]
%   img - The image frame to find peaks in
%   opt - Peak fiding options, as a struct with fields
%     .interpFactor - Factor by which to interpolate frames 
%     .offset       - Morphological offset
%     .threshold    - Minimum normalized peak intensity [0, 1]
%     .sigma        - Gaussian smoothing parameter [pixels]
%     .rgnSize      - Number pixels for valid peak region [min, max]
%     .usePeak      - Whether to use peak region intensity
%                     [false recommended, will use centroid instead]
%
% Outputs
%   xPeaks   - 1-by-N vector of x peak positions [m]
%   zPeaks   - 1-by-N vector of z peak positions [m]
%   peakVals - 1-by-N vector of z peak positions [m]
%
%   
% For further details, please see (and consider citing if used)
%   [1] "Morphological Reconstruction Improves Microvessel Mapping 
%       in Super-Resolution Ultrasound". IEEE T. UFFC 68(6) (2021) 
%       DOI: 10.1109/TUFFC.2021.3057540
%   [2] "MR for ULTRA-SR: Improved Localization with Morphological Image 
%       Processing". IEEE IUS Proc. (2022) 
%       DOI: 10.1109/IUS54386.2022.9957276
%
%         Scott Schoen Jr | MGH-CURT | sschoenjr@mgh.harvard.edu
%
% *************************************************************************

function [xPeaks, zPeaks, peakVals] = isolateBubbles( x, z, img, opt )

% Try to read in from options struct
try
    INTERP_FACTOR = opt.interpFactor;
    MORPH_OFFSET = opt.offset;
    INT_THRESHOLD = opt.threshold; % Minimum normalized peak intensity
    FILTER_SIZE_PX = opt.sigma; % Gaussian filter sigma [px]
    RGN_SIZE = opt.rgnSize; % Num pixels allowed in peak region
    USE_PEAK = opt.usePeak; % Use location of peak intensity instead of centroid
    validOptions = true;
catch
    warning( 'Invalid Options! Using default values.' );
    validOptions = false;
end

% Set parameters
if nargin < 4 || ~validOptions
    INTERP_FACTOR = 2;
    MORPH_OFFSET = 0.25;
    INT_THRESHOLD = 0.5; % Minimum normalized peak intensity
    SIZE_THRESHOLD = 5; % Max region size Pixels
    FILTER_SIZE_PX = 0; % No smoothing
    RGN_SIZE = [0, inf]; % No smoothing
    USE_PEAK = false; % Centroid
end

% Interpolate vectors if desired
if INTERP_FACTOR > 1
    
    % Set flags and store old vectors
    interpolateData = true;
    img0 = img;
    x0 = x;
    z0 = z;
    
    xVec = x0( 1, : );
    dx = xVec(2) - xVec(1);
    zVec = z0( :, 1 );
    dz = zVec(2) - zVec(1);
    
    % Get new position vectors
    xVecInterp = min(x0(:)) : dx./INTERP_FACTOR : max(x0(:));
    dx = xVecInterp(2) - xVecInterp(1);
    zVecInterp = min(z0(:)) : dz./INTERP_FACTOR : max(z0(:));
    dz = zVecInterp(2) - zVecInterp(1);
    [x, z] = meshgrid( xVecInterp, zVecInterp );
    
    newDim = [ size(x,1), size(x,2) ];
    img = imresize( img, newDim, 'cubic' );
    
    xVec = xVecInterp;
    zVec = zVecInterp;
    
else
    x0 = x;
    z0 = z;
    xVec = x0( 1, : );
    dx = xVec(2) - xVec(1);
    zVec = z0( :, 1 );
    dz = zVec(2) - zVec(1);
end

% Smooth image
if FILTER_SIZE_PX > 1
    filtKernel = INTERP_FACTOR.*FILTER_SIZE_PX;
    smoothedImg = imgaussfilt( img, filtKernel );
else
    smoothedImg = img;
end

% Shift image to create marker and mask
mask = smoothedImg;
marker = smoothedImg - MORPH_OFFSET;

% Perform dilation
dilation = imreconstruct(marker, mask);
peaksImg = smoothedImg - dilation;

% --- Detect peaks ---

% First set all values below threshold to 0
pkIm0 = peaksImg;
peaksImg(peaksImg < INT_THRESHOLD.*max(peaksImg(:))) = 0;

% Then binarize the image
bwPeakImg = imbinarize(peaksImg);

% Get properties of non-zero regions
s = regionprops(bwPeakImg, ...
    'PixelIdxList', ...
    'Centroid', ...
    'MajorAxisLength', ...
    'MinorAxisLength', ...
    'Orientation');

% Keep only regions that fit within the size tolerance
sizeRange = INTERP_FACTOR.*RGN_SIZE;
Ns = length(s);
keepInds = zeros(1, Ns );
for sCount = 1 : Ns
    numPixels = length( s(sCount).PixelIdxList );
    bigEnough = numPixels >= sizeRange(1);
    smallEnough = numPixels <= sizeRange(2);
    if bigEnough && smallEnough
        keepInds(sCount) = 1;
    end
end
s0 = s;
s = s(find(keepInds));

% Initialize vectors for this frame
xPeaks = [];
zPeaks = [];
tPeaks = [];
peakVals = [];

% For each region
for peakCount = 1:length(s)   
    
    % Store values for this peak
    if USE_PEAK
        
        % Get pixels for that region
        pixelList = s(peakCount).PixelIdxList;
                
        % Get max intensity
        pixelIntensity = smoothedImg(pixelList);
        [peakValue, maxIdx] = max(pixelIntensity);
        
        % Store
        pixelIdx = pixelList(maxIdx);
        [zIdx, xIdx] = ind2sub(size(img), pixelIdx);
        zPeaks(peakCount) = z(zIdx, xIdx);
        xPeaks(peakCount) = x(zIdx, xIdx);
        peakVals(peakCount) =  peakValue;
        
    else
        % Get centroid
        centroid_inds = round( s(peakCount).Centroid );
        xIdx = centroid_inds(1);
        zIdx = centroid_inds(2);
        zPeaks(peakCount) = z(zIdx, xIdx);
        xPeaks(peakCount) = x(zIdx, xIdx);
        peakVals(peakCount) = smoothedImg(zIdx, xIdx);
    end
    
end

end

