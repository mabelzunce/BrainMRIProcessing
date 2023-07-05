% mode = 0: meanFMRI
% mode = 1:
function SaveOverlayImageWithRois(fmri, atlas, cmap, filename, mode)
    if ndims(fmri) == 4
        meanFmri = mean(fmri,4);
        meanFmri = (meanFmri-min(meanFmri, [], 'all'))./(max(meanFmri, [], 'all')+min(meanFmri, [], 'all'));
    end
    if isempty(cmap)
        cmap = colormap('lines');
    end
    for j = 1 : size(meanFmri, 3)
        % Go through each slice and generate a fused image:
        numLabelsThisSlice = max(max(atlas(:,:,j)));
        imageWithMasks = labeloverlay(1.5*meanFmri(:,:,j)', atlas(:,:,j)', 'Transparency',0.8, 'IncludedLabels',...
            [1:numLabelsThisSlice], 'Colormap',cmap);
        if ndims(imageWithMasks) == 3
            [imind,cm] = rgb2ind(imageWithMasks,256); 
        else
            % Usually means that there are not any label
            imind = imageWithMasks;
            cm = colormap('gray');
        end
        % Write to the GIF File 
        if j == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
        j = j + 1;
    end