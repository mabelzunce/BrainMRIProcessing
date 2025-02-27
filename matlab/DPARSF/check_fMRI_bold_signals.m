% Function that check the fmri signals stored in a matrix.
% roiSignals: matrix that stores in each column the BOLD signal for a ROI.
% (Numbers od rows: length of the BOLD signal. Number of cols: number of
% ROIs.

function fig = check_fMRI_bold_signals(roiSignals) 

    fig = figure('units','normalized','outerposition',[0 0 1 1]) ;
    
    roiSignals(isnan(roiSignals)) = 0;
    % Cross correlation:
    crossCorrSignals = corr(roiSignals) ; % cross corr
    
    % [ SEÑALES BOLD ] 
    subplot(2,3,1) ;
    for j = 1 : size(roiSignals,2)
        plot (roiSignals(:,j)) ;
        title ('Señales BOLD') ;
        xlim( [0 size(roiSignals,1)] ) ;
        hold on ;
    end
    hold off ;
    
    % [ MATRICES ] 
    subplot(2,3,2) ;
    imagesc(crossCorrSignals) ;
    sp7 = subplot(2,3,3) ;
    spy(crossCorrSignals>0.5) ;
    title ('mayor a 0.5') ;
    
    % [ HISTOGRAMA ] 
    subplot(2,3,4) ;
    histogram(crossCorrSignals)
    xlim([-1,1]) ;
    title('histograma correlación cruzada') ;
    
    % [ AUTOCORRELACION ] 
    subplot(2,3,5) ;
    lags= 75 ;
    for j = 1 : size(roiSignals,2)
        y = autocorr(roiSignals(:,j),'NumLags',lags ) ;

        plot (y) ;
        title ('acf') ;
        hold on ;
    end
    hold off ;

    
    % Power spectrum of all signals
    roiSignalsZ = normalize(roiSignals);
    roiSignalsZ(isnan(roiSignalsZ)) = 0;
    % Concatenate:
    roiSingalsZ_concat = reshape(roiSignalsZ, [], 1);

    % Concatenate with interpolation between ROIs.
    % We add a row to the matrix with (signal(roi,end)+signal(roi+1,1))/2
    lastRow = [roiSignalsZ(1, 2:end) roiSignalsZ(1, end)]; % Repeat tthe last column.
    lastRow = (lastRow + roiSignalsZ(end, 1:end))./2;
    roiSignalsInterpZ = [roiSignalsZ; lastRow];
    roiSingalsInterpZ_concat = reshape(roiSignalsInterpZ, [], 1);
    % Show power spectrum:
    subplot(2,3,6) ;
    pspectrum(roiSingalsInterpZ_concat)
end 