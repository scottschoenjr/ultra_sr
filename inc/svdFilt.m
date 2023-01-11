% *************************************************************************
%
% ULTRA-SR | SVD Filtering
%
% Inputs
%   rawData - Stack of acquired US frames [N_frames, Nx, Nz]
%   svdVals - If 2 element vector, specifies minimum and maximum singular
%             values to retain. If a scalar, singular values above this
%             value are kept.
%
% Outputs
%   filteredData - Stack of filtered US frames [N_frames, Nx, Nz]
%   Sraw         - Vector of singular values of raw data
%
% For further details, please see 
%   [1] "Spatiotemporal Clutter Filtering of Ultrafast Ultrasound Data 
%       Highly Increases Doppler and fUltrasound Sensitivity". IEEE T. 
%       Med. Imaging 34(11) (2015)
%       DOI: 10.1109/TMI.2015.2428634
%
%         Scott Schoen Jr | MGH-CURT | sschoenjr@mgh.harvard.edu
%
% *************************************************************************

function [ filteredData, Sraw ] = svdFilt( rawData, svdVals )

[numFrames, numRows, numCols] = size( rawData );
rawData = permute( rawData, [2,3,1] );

% Set default range of SVD values if not specified
if nargin == 2
    if length( svdVals ) == 2
        SVD_VALS = svdVals;
    else
        SVD_VALS = NaN;
        SVD_THRESH = svdVals;
    end
else
    SVD_VALS = NaN;
    SVD_THRESH = 0.1;
end

% Reshape data
A = reshape( rawData, [numRows.*numCols, numFrames] );

% Get decomposition
[ U, S, V ] = svd(A, 'econ');

% Delete largest singular values (i.e., those corresponding to tissue
% motion)
Sraw = S; % To look at later
if ~any( isnan(SVD_VALS) ) && all( SVD_VALS > 0 )
    rmInds = [ 1 : SVD_VALS(1) - 1, SVD_VALS(2) + 1 : numFrames ];
    S( rmInds, rmInds ) = 0;
elseif ~any( isnan(SVD_VALS) ) && all( SVD_VALS < 0 )
    SVD_VALS = -SVD_VALS;
    rmInds = [ SVD_VALS(1) : SVD_VALS(2) ];
    S( rmInds, rmInds ) = 0;
else
    S( S > SVD_THRESH.*max(S(:)) ) = 0;
end

% Reconstruct data
Arecon = U*S*V';

% Reshape back to images
filteredData = reshape( Arecon, [numRows, numCols, numFrames] );
filteredData = permute( filteredData, [3, 1, 2] );

end

