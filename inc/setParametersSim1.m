% Set filtering options depending on simulated dataset 1

% Filtering options
filtOpts.minSVD = 5;
filtOpts.maxSVD = 400;
filtOpts.psfX = 2; % [wvl]
filtOpts.psfZ = 1;
filtOpts.iter = 4;
filtOpts.noise = 0.00;

% Set peak finding options
pkOpts.offset = 0.5;
pkOpts.threshold = 0.95;
pkOpts.interpFactor = 8;
pkOpts.rgnSize = [0.5, 7.5].^2;
pkOpts.sigma = 2;
pkOpts.usePeak = false;
pkOpts.proxThreshold = 0.5; % Must be this close to a true position to be correct [lambda]

% Tracking options
trackOpts.maxV = 25E-3; % [m/s]
trackOpts.minV = 1E-3; % [m/s]
trackOpts.minFrames = 20;

% Image options
imgOpts.interpFactor = 4;
imgOpts.pointDither = 0.25;

% Plotting options
xLimits = 1E-3.*[-12, 12];
yLimits = 1E-3.*[48, 90];

% Image options
imgOpts.interpFactor = 4;
imgOpts.psf_wvl = 1./3;
imgOpts.pointDither = 0.25;
% Insets
imgOpts.inset(1).yLim = 1E-3.*[60, 80];
imgOpts.inset(1).xLim = 1E-3.*[0, 10];

