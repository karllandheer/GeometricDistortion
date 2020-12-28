%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
    function [projAbsCorr,projFieldCorr,f_succ] = FAMASITO_GeoCorrection(projAbsOrig,projFieldOrig,fBackgr2,varargin)
%% 
%%  FASTMAP Shim Tool (FAMASITO) - Spatial distortion correction.
%%  Correction of spatial distortion induced by B0 background field along
%%  readout direction in projection-based FASTMAP-type B0 field mapping.
%%
%%
%% --- data despription ---
%% nspecC:  number of points along projection readout
%% nProj:   number of projections (here: 6)
%% nTE:     number of echo time delays (here: 2)
%%
%%
%% --- function inputs ---
%% projAbsOrig: measured projections potentially including spatial distortions
%% format: nspecC, nProj, nTE
%%
%% projFieldOrig: measured field projections potentially including spatial distortions
%% format: nspecC, nProj
%%
%% fBackgr2: projection-specfic background field distribution. This is the
%% source of spatial distortion. To this end, the sum of the readout gradient
%% and this background field distribution is mapped back onto the gradient 
%% itself and then employed for spatial correction of the projection
%% geomtries. Note that the background field distribution is provided at
%% twice the FOV in order to avoid boarder effects.
%% format: nspecC, nProj
%%
%% projAbsRef: reference projection for spatial shift correction (as opposed
%% to a stretch correction). Typical example: offset from Z2. A spatial 
%% shift correction is applied if a 4th function
%% argument is assigned.
%%
%%
%% --- function output ---
%% projAbsCorr: projections after spatial corrections
%% format: nspecC, nProj, nTE
%%
%% projFieldCorr: field projections after spatial corrections
%% format: nspecC, nProj
%%
%% f_succ: success flag 
%% format: binary (0: failed, 1: success)
%% 
%%
%%  Karl Landheer & Christoph Juchem
%%  10/2017
%%  Columbia University
%%  All rights reserved.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global shim flag

FCTNAME = 'FP_Shim_GeoCorrection';


%--- init success flag ---
f_succ = 0;
            
%--- argument handling ---
if nargin==3            % stretch correction only
    f_shiftCorr  = 0;
elseif nargin==4        % stretch and shift correction
    f_shiftCorr  = 1;
    projShiftRef = FP_Check4NumR( varargin{1} );
else
    fprintf('%s ->\n%i function arguments are not supported. Program aborted.\n',FCTNAME,nargin)
    return
end


%--- double FOV handling for spatial remapping ---
% note the difference to shim.posVecCm2 to assure identical frequency
% positions at the center 1x FOV piece.
HzPPt        = shim.sw_h/(shim.nspecC-1);                       % Hz per point in original FOV = bandwidth (to be applied to both 1xFOV and 2xFOV)
limHz2       = shim.sw_h/2 + shim.nspecC/2*HzPPt;               % limHz2: frequency limits of 2xFOV bandwidth
fVecTarget2  = -limHz2:2*limHz2/(2*shim.nspecC-1):limHz2;       % 2x readout gradient (not used)
% here the point shim.nspecC+1 of the 2xFOV bandwidth is identical to
% the frequency of the first point of the original 1xFOV bandwidth
ind2First    = round(shim.nspecC/2)+1;                          % first index of center 1xFOV piece in larger 2xFOV vector
ind2Last     = ind2First+shim.nspecC-1;                         % last index of center 1xFOV piece in larger 2xFOV vector
fVecReal2    = zeros(2*shim.nspecC,shim.nProj);                 % init real apparent field matrix as target + background

%--- magnitude projection handling ---
projAbsOrig2 = zeros(2*shim.nspecC,shim.nProj,2);               % init original projection matrix at twice the FOV
projAbsOrig2(ind2First:ind2Last,:,:) = projAbsOrig;             % assign to center position
projAbsCorr2 = zeros(2*shim.nspecC,shim.nProj,2);               % init correction projection matrix

%--- field projection handling ---
projFieldOrig2 = zeros(2*shim.nspecC,shim.nProj);               % init original projection matrix at twice the FOV
projFieldOrig2(ind2First:ind2Last,:) = projFieldOrig;           % assign to center position
projFieldCorr2 = zeros(2*shim.nspecC,shim.nProj);               % init correction projection matrix

%--- stretch correction ---
for pCnt = 1:shim.nProj
    fVecReal2(:,pCnt)      = fVecTarget2 + fBackgr2(:,pCnt).'/2;    % factual frequency distribution
    projAbsCorr2(:,pCnt,1) = pchip(fVecReal2(:,pCnt),projAbsOrig2(:,pCnt,1),fVecTarget2);
    projAbsCorr2(:,pCnt,2) = pchip(fVecReal2(:,pCnt),projAbsOrig2(:,pCnt,2),fVecTarget2);
    projFieldCorr2(:,pCnt) = pchip(fVecReal2(:,pCnt),projFieldOrig2(:,pCnt),fVecTarget2);
end

%--- shift correction ---
if f_shiftCorr
    % The shift correction is determined on the original FOV projections 
    % and then applied to the extended FOV data before the original FOV is
    % extracted. This is possible as they both use the identical spatial
    % resolution. It is based on the assumption, however, that signals 
    % remain well-defined at the FOV borders, i.e. the object extends 
    % beyond the FOV.  
    
    %--- determine optimized shift correction ---
    % note:
    % 1) the first TE magnitude projections after unstretching are used
    % 2) the resultant resolution of this approach is limited to the point resolution
    optMeasure = zeros(shim.nspecC,shim.nProj);            
    for nCnt = 1:shim.nspecC
        optMeasure(nCnt,:) = sum(abs(circshift(projAbsCorr2(ind2First:ind2Last,:,1),[nCnt 0 0]) - projShiftRef));
    end
    [optMeasVal,optShift] = min(optMeasure);
    % center
    leftInd = find(optShift>shim.nspecC/2);
    optShift(leftInd) = optShift(leftInd)-shim.nspecC;
    
    %--- apply spatial shift correction to projections ---
    projAbsOrigShift   = zeros(shim.nspecC,shim.nProj,shim.nTE);        % init
    projFieldOrigShift = zeros(shim.nspecC,shim.nProj,shim.nTE);        % init
    for pCnt = 1:shim.nProj
        % single FOV
        projAbsOrigShift(:,pCnt,1) = circshift(projAbsOrig(:,pCnt,1),optShift(pCnt));
        projAbsOrigShift(:,pCnt,2) = circshift(projAbsOrig(:,pCnt,2),optShift(pCnt));
        projFieldOrigShift(:,pCnt) = circshift(projFieldOrig(:,pCnt),optShift(pCnt));
    
        % double FOV
        projAbsCorr2(:,pCnt,1) = circshift(projAbsCorr2(:,pCnt,1),optShift(pCnt));
        projAbsCorr2(:,pCnt,2) = circshift(projAbsCorr2(:,pCnt,2),optShift(pCnt));
        projFieldCorr2(:,pCnt) = circshift(projFieldCorr2(:,pCnt),optShift(pCnt));
    end
    
    %--- result visualization ---
    % shown is the 1xFOV data since that is used to determine the shifts
    if flag.verbose
        fhShiftAna = figure;
        set(fhShiftAna,'NumberTitle','off','Color',[1 1 1],'Position',[110 120 1700 875],...
            'Name',sprintf('FAMASITO Shift Correction'),'Tag','Calib')

        for pCnt = 1:shim.nProj
            %--- optimization shift curves ---
            subplot(2,shim.nProj,pCnt)
            plot(optMeasure(:,pCnt))
            [plotLim(1) plotLim(2) plotLim(3) plotLim(4)] = FP_IdealAxisValues(optMeasure(:,pCnt));    
            axis(plotLim)
            if pCnt==1
                ylabel('optimization measure')
            end
            switch pCnt
                case 1
                    title(sprintf('Trace #%.0f: +X/+Y',pCnt))
                case 2
                    title(sprintf('Trace #%.0f: -X/+Y',pCnt))
                case 3
                    title(sprintf('Trace #%.0f: +X/+Z',pCnt))
                case 4
                    title(sprintf('Trace #%.0f: -X/+Z',pCnt))
                case 5
                    title(sprintf('Trace #%.0f: +Y/+Z',pCnt))
                case 6
                    title(sprintf('Trace #%.0f: -Y/+Z',pCnt))
                    xlabel('position [cm]')
            end

            %--- magnitude projections ---
            subplot(2,shim.nProj,shim.nProj+pCnt)
            hold on
            plot(shim.posVecCm,projAbsOrig(:,pCnt,1),'r')
            plot(shim.posVecCm,projShiftRef(:,pCnt))
            plot(shim.posVecCm,projAbsOrigShift(:,pCnt),'g')
            hold off
            if pCnt==shim.nProj
                legend('Orig','Ref','Corr')
            end
            [plotLim(1) plotLim(2) plotLim(3) plotLim(4)] = FP_IdealAxisValues(shim.posVecCm,projAbsOrig(:,pCnt,1));    
            axis(plotLim)
            if pCnt==1
                ylabel('amplitude [a.u.]')
            end
        end
    end
end

%--- regular FOV handling: extraction of original FOV ---
% as center pieces of extended vectors
fVecTarget    = fVecTarget2(ind2First:ind2Last);                % readout gradient
fVecReal      = fVecReal2(ind2First:ind2Last,:);
projAbsCorr   = projAbsCorr2(ind2First:ind2Last,:,:);
projFieldCorr = projFieldCorr2(ind2First:ind2Last,:);

%--- verbose printout ---
if flag.verbose
    fhCorrAna = figure;
    set(fhCorrAna,'NumberTitle','off','Color',[1 1 1],'Position',[110 120 1700 875],...
        'Name',sprintf('FAMASITO Individual Geometry Correction'),'Tag','Shim')

    for pCnt = 1:shim.nProj
        %--- field traces ---
        subplot(3,shim.nProj,pCnt)
        hold on
        plot(shim.posVecCm,fVecReal(:,pCnt),'r')
        plot(shim.posVecCm,fVecTarget)
        if pCnt==shim.nProj
            legend('Orig','Corr')
        end
        [minX maxX minY maxY] = FP_IdealAxisValues(shim.posVecCm,fVecReal(:,pCnt));
        if flag.shimPlotRange==1            % full
            axis([minX maxX minY maxY])
        elseif flag.shimPlotRange==2        % ROI
            [minXroi maxXroi minYroi maxYroi] = ...
                FP_IdealAxisValues(shim.posVecCm(shim.roiLimInd(pCnt,1):shim.roiLimInd(pCnt,2)),...
                                   fVecReal(shim.roiLimInd(pCnt,1):shim.roiLimInd(pCnt,2),pCnt));
            rgExt = 0.05*(maxYroi-minYroi);        % 5% extension
            axis([minX maxX minYroi-rgExt maxYroi+rgExt])
        else                            % selected range
            axis([minX maxX shim.plotRgManMin shim.plotRgManMax])
        end        
                
        %--- visualization of voxel limits ---
        if flag.shimPlotVoxel
            if flag.shimPlotRange==1            % full
                plot([shim.roiLimCm(pCnt,1) shim.roiLimCm(pCnt,1)],[min(minY) max(maxY)],'g')
                plot([shim.roiLimCm(pCnt,2) shim.roiLimCm(pCnt,2)],[min(minY) max(maxY)],'g')
            elseif flag.shimPlotRange==2        % ROI
                plot([shim.roiLimCm(pCnt,1) shim.roiLimCm(pCnt,1)],[minYroi-rgExt maxYroi+rgExt],'g')
                plot([shim.roiLimCm(pCnt,2) shim.roiLimCm(pCnt,2)],[minYroi-rgExt maxYroi+rgExt],'g')
            else                                % selected range
                plot([shim.roiLimCm(pCnt,1) shim.roiLimCm(pCnt,1)],[shim.plotRgManMin shim.plotRgManMax],'g')
                plot([shim.roiLimCm(pCnt,2) shim.roiLimCm(pCnt,2)],[shim.plotRgManMin shim.plotRgManMax],'g')
            end
        end
        hold off
        
        if pCnt==1
            ylabel('readout field [Hz]')
        end
        switch pCnt
            case 1
                title(sprintf('Trace #%.0f: +X/+Y',pCnt))
            case 2
                title(sprintf('Trace #%.0f: -X/+Y',pCnt))
            case 3
                title(sprintf('Trace #%.0f: +X/+Z',pCnt))
            case 4
                title(sprintf('Trace #%.0f: -X/+Z',pCnt))
            case 5
                title(sprintf('Trace #%.0f: +Y/+Z',pCnt))
            case 6
                title(sprintf('Trace #%.0f: -Y/+Z',pCnt))
        end
        
        %--- magnitude traces ---
        subplot(3,shim.nProj,shim.nProj+pCnt)
        hold on
        plot(shim.posVecCm,projAbsOrig(:,pCnt,1),'r')
        plot(shim.posVecCm,projAbsCorr(:,pCnt,1))
        hold off
        if pCnt==shim.nProj
            legend('Orig','Corr')
        end
        [plotLim(1) plotLim(2) plotLim(3) plotLim(4)] = FP_IdealAxisValues(shim.posVecCm,projAbsOrig(:,pCnt,1));
        axis(plotLim)
        if pCnt==1
            ylabel('amplitude [a.u.]')
        end
        
        %--- field traces ---
        subplot(3,shim.nProj,2*shim.nProj+pCnt)
        hold on
        plot(shim.posVecCm,projFieldOrig(:,pCnt),'r')
        plot(shim.posVecCm,projFieldCorr(:,pCnt))
        hold off
        if pCnt==shim.nProj
            legend('Orig','Corr')
        end
        [plotLim(1) plotLim(2) plotLim(3) plotLim(4)] = FP_IdealAxisValues(shim.posVecCm,projFieldOrig(:,pCnt));
        axis(plotLim)
        xlabel('position [cm]')
        if pCnt==1
            ylabel('background field [Hz]')
        end
    end
end

%--- update success flag ---
f_succ = 1;






