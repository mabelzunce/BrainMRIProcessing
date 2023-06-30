function fig = check_fMRI_DPARSF (working_dir ,  T , subID , signals_dir , outputPath, save ) 
    
    % VARIABLES:
    % carpeta : (char) nombre de la carpeta de DPARFS con las nifti corregidas
    % ej carpeta = 'FunImgARWSDC/' 
    % working_dir : (char)  es el de DPARFS 
    % signals_dir : (char) es el Results/nombre de la carpeta con las señales 
    % ej 'Results/ROISignals_FunImgARWSDCF/'
    % T : tabla con los datos del sujeto (table) 
    % subID (char)    
    % save = 1 guarda las figuras 
    % save = 0 no guarda las figuras 
    
    
    figure('Name',  subID ,'units','normalized','outerposition',[0 0 1 1]) ;   
    fig = gcf;
    
    
    % [ IMAGENES ] 
    p_pictures = strcat(working_dir,'PicturesForChkNormalization/',subID,'.tif');
    pic = imread(p_pictures) ;
    sp1 = subplot(3,3,1) ;
    imshow(pic) ;   
     
    % [ DATA ] 
    indice = string(T.SubjectID) == string(subID) ; 
    indice = find(indice,1) ;
    escaner = 'Philips' ; % usar  strsplit con delimiter ;
    estudio = T.Phase(indice) ; estudio = char(estudio) ;
    genero = T.Sex(indice) ; genero = char(genero) ;
    edad = T.Age(indice) ; edad = string(edad) ;
    fecha = T.StudyDate(indice) ; fecha = char(fecha) ;
    grupo = T.ResearchGroup(indice) ; grupo = char(grupo) ;
    %TP = T.TimePoints(indice)-5 ; TP = string(TP) ;
    %procesado = T.Procesado(indice) ; procesado = string(procesado) ;

    % [ CORRELACIONES ] 
    ROI = strcat( 'ROISignals_' , subID , '.mat' ) ;
    p_ROI = strcat(working_dir, signals_dir , ROI) ;
    S = load (p_ROI) ;
    S = S.ROISignals ; % señales BOLD   
    corr_S = corr(S) ; % cross corr
    
    % [ SEÑALES BOLD ] 
    for j = 1:116
        sp2 = subplot(3,3,2) ;
        plot (S(:,j)) ;
        title ('Señales BOLD') ;
        xlim( [0 size(S,1)] ) ;
        hold on ;
    end
    hold off ;
    
    % [ MATRICES ] 
    sp4 = subplot(3,3,4) ;
    imagesc(corr_S) ;
    sp7 = subplot(3,3,7) ;
    spy(corr_S>0.5) ;
    title ('mayor a 0.5') ;
    
    % [ HISTOGRAMA ] 
    sp3 = subplot(3,3,3) ;
    histogram(corr_S)
    xlim([-1,1]) ;
    title('histograma correlación cruzada') ;
    
    % [ AUTOCORRELACION ] 
    lags= 75 ;
    for j = 1:116
        y = autocorr ( S(:,j),'NumLags',lags ) ;
        sp5 = subplot(3,3,5) ;
        plot (y) ;
        title ('acf') ;
        hold on ;
    end
    hold off ;
    
%     % [ MEDIA VS VARIANZA ] 
%     m=[] ; v=[] ;
%     for j = 1:116
%         m(j) = abs(mean ( S(:,j) )) ;
%         v(j) = abs(var ( S(:,j) )) ;
%     end
%     sp6 = subplot(3,3,6) ;
%     loglog(m,v,'r*') ;
%     xlabel('media') ;
%     ylabel('varianza') ;   
    
    roiSignalsZ = normalize(S);
    % Concatenate:
    roiSingalsZ_concat = reshape(roiSignalsZ, [], 1);
    % Show power spectrum:
    subplot(3,3,6) ;
    pspectrum(roiSingalsZ_concat)

    % [ MOVIMIENTO ]
    report_files = dir ( strcat(working_dir,'/RealignParameter/' ,subID)) ;
    for i=3:length(report_files) % busca el archivo
        if contains( string(report_files(i).name) , 'rp')
            nombre = string(report_files(i).name)  ;
            report = load(strcat(working_dir,'/RealignParameter/' ,subID,'/', nombre)) ;
        end
    end
    sp8 = subplot(3,3,8) ;
    plot(report(:,1)) ; 
    hold on ; 
    plot(report(:,2)) ; 
    plot(report(:,3)) ;
    title('movimiento') ;
    legend({'x','y','z'},'Location','southeast','FontSize',6) ;
    xlim([ 0 , size(report,1)] ) ;
    
    
    % [ DATOS DEL SUJETO ]
    sp9 = subplot(3,3,9) ;
    %title('Datos del sujeto') ;
    subID = string(subID) ;
    text(0.1,0.7,subID,'FontSize',12) ;
    text(0.55,0.7,grupo,'FontSize',12) ;
    text(0.7,0.7,estudio,'FontSize',12) ;
    text(0.1,0.35,strcat('edad: ', edad),'FontSize',12);
    text(0.55,0.35,genero,'FontSize',12) ;
    text(0.7,0.35,fecha,'FontSize',12) ;
    %text(0.1,0.0,strcat('procesado:',procesado') ,'FontSize',12);
    %text(0.55,0.0,strcat('Time Points =',TP),'FontSize',12) ;
    axis off ;
    
    
    % [ POSICIONES ] 
    % [ left , bottom , width , heigth ]
    % Distance from the left edge 
    % Distance from the bottom edge
    % Distance between the right and left inner edges of the figure
    % Distance between the top and bottom inner edges of the window.
    
%     mgleft = -0.1 ; mvrigth = 0.075 ;
%     mgbottom = -0.05 ;
%     wdleft = 0 ; wdcenter = 0.085 ; wdrigth = -0.01 ;
%     
%     
%     pos1 = get(sp1, 'Position') ;
%     new_pos1 = pos1 + [0 0 0 0] 
%     set(sp1, 'Position',new_pos1 ) ; 
%     
%     pos2 = get(sp2, 'Position') ;
%     new_pos2 = pos2 + [0 0 wdcenter 0] 
%     set(sp2, 'Position',new_pos2 ) ; 
%     
%     pos3 = get(sp3, 'Position') ;
%     new_pos3 = pos3 + [mvrigth 0 wdrigth 0] 
%     set(sp3, 'Position',new_pos3 ) ; 
%     
%     pos4 = get(sp4, 'Position') ;
%     new_pos4 = pos4 + [mgleft 0 wdleft 0] 
%     set(sp4, 'Position',new_pos4 ) ; 
%     
%     pos5 = get(sp5, 'Position') ;
%     new_pos5 = pos5 + [0 0 wdcenter 0] 
%     set(sp5, 'Position',new_pos5 ) ;
%     
%     pos6 = get(sp6, 'Position') ;
%     new_pos6 = pos6 + [mvrigth 0 wdrigth 0] 
%     set(sp6, 'Position',new_pos6 ) ;
%     
%     pos7 = get(sp7, 'Position') ;
%     new_pos7 = pos7 + [mgleft mgbottom wdleft 0] 
%     set(sp7, 'Position',new_pos7 ) ;
%     
%     pos8 = get(sp8, 'Position') ;
%     new_pos8 = pos8 +[0 mgbottom wdcenter 0] 
%     set(sp8, 'Position',new_pos8 ) ;
%     
%     pos9 = get(sp9, 'Position') ;
%     new_pos9 = pos9 +[mvrigth-0.03 mgbottom wdrigth 0] 
%     set(sp9, 'Position',new_pos9 ) ;
    
    % [ GUARDAR ARCHIVOS ] 
    if save==1
        if ~exist(outputPath, 'dir')
           mkdir(outputPath)
        end
        nombre = strcat(subID) ;
        %savefig(nombre)
        saveas(fig,[outputPath char(subID) ],'epsc')
    end
    
    
   fig = gcf;
   
end 