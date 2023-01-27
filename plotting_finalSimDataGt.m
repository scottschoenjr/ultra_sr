% Plotting Script - Run after processing

close all;
clc;

PLOT_TRACKING = false;

% Set inset number
insetNum = 1;
dataFraction = [0.00, 0.02];

% Set frame Nos. to Plot (if the whole dataset was prcessed)
if ( DATA_FRACTION(2) - DATA_FRACTION(1) ) >= 1
        
    % Save full sets
    rawPeaks = peakMatMR;
    rawTracks = trackStructMr;
    rawFrames = numFrames;
    
    % Get frames of time range of interest
    Nf = size( bwVideo, 1 );
    startFrame = round( dataFraction(1).*Nf + 1);
    endFrame = round( dataFraction(2).*Nf - 1);
    numFrames = endFrame - startFrame + 1;
    
    % Keep only those in desired region
    frameVals = peakMatMR( :, 1 );
    keepInds = frameVals >= startFrame & frameVals <= endFrame;
    peakMatMR = rawPeaks( keepInds, :, : );
    trackStructMr = rawTracks( startFrame : endFrame );

else
    
    % Otherwise just keep it all
    numFrames = size(bwVideo, 1);
    
end



% Retain peaks only in certain range of amplitudes
peakAmps = peakMatMR( :, 6 );
minAmp = mean(peakAmps) - std(peakAmps);
maxAmp = mean(peakAmps) + 1E6.*std(peakAmps);

keepInds = peakAmps >= minAmp & peakAmps <= maxAmp;
peakWeights = peakAmps./max(abs(peakAmps(:)));

%% Make Image
figure();

% Interpolate image to this grid
xSr = imresize(x, imgOpts.interpFactor);
zSr = imresize(z, imgOpts.interpFactor);

% Set PSF Size
psfSize = imgOpts.psf_wvl.*lambda;

[img] = makeSrImg( peakMatMR(keepInds, 2), peakMatMR(keepInds, 3), ...
    psfSize, xSr, zSr, false, peakWeights(keepInds) );
imPlot = 20.*log10( abs(img)./max(abs(img(:))) );

% Plot full image
pcolor( 1E3.*xSr, 1E3.*zSr, imPlot );

colormap hot;
shading flat;
caxis([-32, 0]);

axis equal;
axis tight;
axis ij;

xlabel( '[mm]' );
xlim(1E3.*xLimits);
ylim(1E3.*yLimits);

% And inset
figure();
pcolor( 1E3.*xSr, 1E3.*zSr, imPlot );

colormap hot;
shading flat;
caxis([-32, 0]);

axis equal;
axis tight;
axis ij;

xlim(1E3.*imgOpts.inset(insetNum).xLim);
ylim(1E3.*imgOpts.inset(insetNum).yLim);

xlabel( '[mm]' );

if ~PLOT_TRACKING
    return;
end

%% Tracking

maxExpectedVelocity = trackOpts.maxV; % [m/s]
maxDistance = maxExpectedVelocity.*(1./vid.FrameRate); % [m]

% Get cell array of tracks
[tracks, adjTracks, A] = ...
    simpletracker(trackStructMr, 'MaxLinkingDistance', maxDistance );

%% Plot velocities

DeltaT = 1./( ImSet.frame_rate_Hz );

% Initialize
xq = [];
zq = [];
vxq = [];
vzq = [];

for trackCount = 1 : length( tracks )
    
    % Make sure track is long enough
    loopTrack = tracks{trackCount};
    pathTooShort = length( find(~isnan(loopTrack)) ) < trackOpts.minFrames;
    if pathTooShort
        continue;
    end
    
    % Now get position and velocity at each point
    xLoop = [];
    zLoop = [];
    tLoop = [];
    for frameCount = 1 : numFrames
        ld = trackStructMr{frameCount};
        if isnan( loopTrack(frameCount) )
            continue;
        end
        xLoop = [ xLoop, ld( loopTrack(frameCount), 1 ) ];
        zLoop = [ zLoop, ld( loopTrack(frameCount), 2 ) ];
        tLoop = [ tLoop, (frameCount-1).*DeltaT ];
    end
    
    % Get velocity at each point
    vxLoop = diff(xLoop)./diff(tLoop);
    vzLoop = diff(zLoop)./diff(tLoop);
    
    % Store to quiver plot
    xq = [xq, xLoop];
    zq = [zq, zLoop];
    vxq = [vxq, vxLoop(1), vxLoop];
    vzq = [vzq, vzLoop(1), vzLoop];
    
end

% Remove outliers
nSigma = 2;
vMag = sqrt( vxq.^(2) + vzq.^(2) );
vTheta = atan2( vzq, vxq );
meanV = mean( vMag );
stdV = std( vMag );

rmInds = ...
    vMag > trackOpts.maxV | ...
    vMag < trackOpts.minV | ...
    vMag > meanV + nSigma.*stdV;

xq(rmInds) = [];
zq(rmInds) = [];
vxq(rmInds) = [];
vzq(rmInds) = [];
vMag(rmInds) = [];
vTheta(rmInds) = [];

% Get velocity plot
Fv = scatteredInterpolant( xq', zq', vMag' );
weight = ( img./max(img(:)) ).^(1.5);
vMap = weight.*Fv( xSr, zSr );

% Plot inerpolated velocity field
figure();

vmp = vMap;
vmp( vmp < 0.05E-3 ) = NaN;
pcolor( 1E3.*xSr, 1E3.*zSr, 1E3.*vmp );
shading interp;
colormap jet;

axis ij;
axis equal;
axis tight;

xlim(1E3.*xLimits);
ylim(1E3.*yLimits);
caxis([0, 1]); % [mm/s]
cbh = colorbar;
ylabel( cbh, 'Speed [mm/s]' );

% Get velocity angle plot
Fa = scatteredInterpolant( xq', zq', vTheta' );
vMap = Fa( x, z );

% Upscale
xVa = imresize(x, 4);
zVa = imresize(z, 4);
vMap = imresize(vMap, 4);
vWeightA = imresize( weight, 4 );

vMap( weight < 0.01 ) = NaN;

figure();
pcolor( 1E3.*xVa, 1E3.*zVa, phasewrap(vMap - pi./2) );
phasemap(16, 'rad' );
shading interp;

axis ij;
axis equal;
axis tight;

xlim(1E3.*xLimits);
ylim(1E3.*yLimits);

%% Color plot

figure();
axes();
ax2 = gca();
set( gca, 'ColorMap', flipud(redblue(16)) );
quiverC2D( 1E3.*xq', 1E3.*zq', 1E3.*vxq', 1E3.*vzq', ...
    'LineWidth', 2.5, ...
    'MaxHeadSize', 0.25, ...
    'DopplerStyle', true, ...
    'ColorBar', true, ...
    'ColorLims', [-10, 10], ...
    'Scale', 0.2 );

axis ij;
axis equal;
axis tight;

xlim(1E3.*xLimits);
ylim(1E3.*yLimits);

%% Colored Arrows

figure();
axes();

% Set limits first for this this color quiver function
axis equal;
axis tight;
axis ij;

xlim(1E3.*imgOpts.inset(insetNum).xLim);
ylim(1E3.*imgOpts.inset(insetNum).yLim);


FUN_quiver_by_plotV2_cmap_patch( 1E3.*xq, 1E3.*zq, 1E3.*vxq, 1E3.*vzq, ...
    0.02, ...
    'zval', phasewrap(-vTheta), ...
    'fill_head', true, ...
    'head_length', 15, ...
    'LineWidth', 1);

cm = phasemap( 256, 'rad' );
Nc = size( cm, 1 );
cm_hsv = rgb2hsv( cm );
cm_hsv( :, 2 ) = min( [1.5.*cm_hsv(:, 2), ones(Nc,1)], [], 2 );
cm_hsv( :, 3 ) = min( [1.2.*cm_hsv(:, 3), ones(Nc,1)], [], 2 );
cm_rgb = hsv2rgb( cm_hsv );
colormap( cm_rgb );

phasecirc( 'location', 'sw' );

%% Plot tracks directly

numTracks = numel(tracks);
trackColors = rand(numTracks, 3);
trackColors = hot(numTracks);

allTrackPoints = vertcat(trackStructMr{:});

figure();
hold all;

for trackCount = 1 : numTracks

    % We use the adjacency tracks to retrieve the points coordinates. It
    % saves us a loop.

    loopTrack = adjTracks{trackCount};
    loopTrackPoints = allTrackPoints(loopTrack, :);

    plot( ...
        1E3.*loopTrackPoints(:,1), ...
        1E3.*loopTrackPoints(:,2), ...
        'Color', 'k', ...
        'LineWidth', 2);

end

axis equal;
axis tight;
axis ij;

xlim(1E3.*imgOpts.inset(insetNum).xLim);
ylim(1E3.*imgOpts.inset(insetNum).yLim);

%% Re-assign
peakMatMR = rawPeaks;
trackStructMr = rawTracks;
numFrames = rawFrames;

%% Display Stats
clc;

fprintf( '---------- MR Statistics (h = %01.2f) ----------\n', ...
    pkOpts.offset );
fprintf( 'Peaks Per Frame: %2.1f +/- %2.1f\n', ...
    mean([frameStats(:).numPeaksMr]), std([frameStats(:).numPeaksMr]) );
fprintf( '    Comp. Time: %3.1f +/- %2.1f ms\n\n', ...
    1E3.*mean([frameStats(:).peakFindingTime]), ...
    1E3.*std([frameStats(:).peakFindingTime]) );

