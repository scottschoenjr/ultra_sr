% Set filtering options depending on experimental dataset

% Filtering options
filtOpts.minSVD = 4;
filtOpts.maxSVD = 400;
filtOpts.psfX = 2; % [wvl]
filtOpts.psfZ = 1;
filtOpts.iter = 4;
filtOpts.noise = 0.00;

% Set peak finding options
pkOpts.offset = 0.4;
pkOpts.threshold = 0.95;
pkOpts.interpFactor = 8;
pkOpts.rgnSize = [0.01, 100].^2;
pkOpts.sigma = 2;
pkOpts.usePeak = false;

% Tracking options
trackOpts.maxV = 25E-3; % [m/s]
trackOpts.minV = 1E-3; % [m/s]
trackOpts.minFrames = 20;

% Plotting options
xLimits = 1E-3.*[-10, 10];
yLimits = 1E-3.*[10, 30];

% Image options
imgOpts.interpFactor = 4;
imgOpts.psf_wvl = 1./3;
imgOpts.pointDither = 0.25;

% Insets
imgOpts.inset(1).xLim = 1E-3.*[-3, 3]; % 0 to 0.1
imgOpts.inset(1).yLim = 1E-3.*[12, 18];

imgOpts.inset(2).xLim = 1E-3.*[-6, -3]; % 0.1 to 0.3
imgOpts.inset(2).yLim = 1E-3.*[18, 21];