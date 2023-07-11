function fig = check_fMRI_DPARSF (working_dir , subID , signals_dir , output, save , T, MERGE ,imgQ_ADNI2,imgQ_ADNI3) 
    
    
    % VARIABLES:
    % carpeta : (char) nombre de la carpeta de DPARFS con las nifti corregidas
    % ej carpeta = 'FunImgARWSDC/' 
    % working_dir : (char)  es el de DPARFS 
    % signals_dir : (char) es el Results/nombre de la carpeta con las señales 
    % ej 'Results/ROISignals_FunImgARWSDCF/'
    % T : tabla con los datos del sujeto (table) 
    % MERGE : tabla de ADNIMERGE (table) 
    % las otras tablas son las de image quality 
    % subID (char)    
    % save = 1 guarda las figuras 
    % save = 0 no guarda las figuras 
    
    
    figure('Name',  subID ,'units','normalized','outerposition',[0 0 1 1]) ;   
    set(gcf,'color','white')
    fig = gcf;
    
    % [ TITULO ]

    %text(0 , 0 , 'text','FontSize',24,'Color','Red') ; 
    
    % [ IMAGENES ] 
    p_pictures = strcat(working_dir,'PicturesForChkNormalization/',subID,'.tif');
    pic = imread(p_pictures) ;
    sp1 = subplot(2,5,1) ;
    imshow(pic) ;  
    hold off ;
    
    % [ DATA ] 
    
    % --> T
    indice = string(T.SubjectID) == string(subID) ; 
    indice = find(indice,1) ;
    imageID = T.ImageID(indice) ; 
    escaner = 'Philips' ; % usar  strsplit con delimiter ;
    estudio = T.Phase(indice) ; estudio = char(estudio) ;
    genero = T.Sex(indice) ; genero = char(genero) ;
    edad = T.Age(indice) ; edad = string(edad) ;
    fecha = T.StudyDate(indice) ; 
    grupo = T.ResearchGroup(indice) ; grupo = char(grupo) ;
    TP = T.TimePoints(indice)-5 ; TP = string(TP) ;
    
    
    % --> image Quality
    
    inicio_ADNI3 = datetime('01/01/2017','InputFormat','MM/dd/uuuu') ;
    fecha = datetime(fecha,'InputFormat','MM/dd/yyyy') ;
    
    if fecha > inicio_ADNI3
        % uso la tabla de ADNI3
        fila = string(imgQ_ADNI3.LONI_IMAGE) == string(imageID) ;
        imgQ = imgQ_ADNI3.SERIES_QUALITY(fila) ;
        study_com = string( imgQ_ADNI3.STUDY_COMMENTS(fila) ) ;
        
    else
        % uso la tabla de ADN2
        fila = string(imgQ_ADNI2.loni_image) == string(imageID) ;
        imgQ = imgQ_ADNI2.series_quality(fila) ;
        study_com = string( imgQ_ADNI2.study_comments(fila) );
    end
    
    
    % --> tests cognitivos:
    usID = str2double ( string(subID(7:10)) ) ;
    img_date = datetime(fecha,'InputFormat','MM/dd/yyyy') ;
    
    pos = find( MERGE.RID==usID ) ;
    % el vector pos contiene el numero de fla (indice) 
    % en el que aparece el sujeto (puede haber mas de uno)
    % por ej, el sujeto 130_S_4641 aparece en las filas 
    % 4450 , 4451, 4452, 4453, 4454, y 4455
    % esos son los numeros que contiene el vector pos

    if ~isempty(pos)   % verifica que el vector no este vacio            
        
        % empieza asumiendo que los datos correctos estan en
        % la primera vez que aparece el sujeto en la tabla, en pos(1)
        % (por ej, el sujeto 130_S_4641 aparece por primera vez en 4450)
        ind = pos(1) ; 
        
        % mindiff es la diferencia de fecha entre la resonancia y los tests
        mindif = abs( caldays( between( img_date , MERGE.EXAMDATE(pos(1)) , 'days')  ));
        tieneMMSE = ~ isnan( MERGE.MMSE(ind) ) ;
        tieneMOCA = ~ isnan( MERGE.MOCA(ind) ) ;
        if tieneMMSE
            MMSE = MERGE.MMSE(ind) ;
        else 
            MMSE = NaN ;
        end
        if tieneMOCA
            MOCA = MERGE.MOCA(ind) ;
        else 
            MOCA = NaN ;
        end
        
        % compara con los otros datos del sujeto que hay en la tabla       
        for j=2:length(pos)
            exam_date = MERGE.EXAMDATE(pos(j)) ;
            ind = pos(j) ;
            dif = abs( caldays( between(img_date,exam_date,'days') ));
            
            % si encuentra una fecha mas cercana a la resonancia
            % toma esos datos
            if dif < mindif
                mindif = dif   ;
                ind = pos(j) ;
                
                tieneMMSE = ~ isnan( MERGE.MMSE(ind) ) ;
                tieneMOCA = ~ isnan( MERGE.MOCA(ind) ) ;
                if tieneMMSE
                    MMSE = MERGE.MMSE(ind) ;
                else 
                    MMSE = NaN ;
                end
                if tieneMOCA
                    MOCA = MERGE.MOCA(ind) ;
                else 
                    MOCA = NaN ;
                end
                
            end
            
        end
    
    else % isempty(pos)=1
        MMSE = 'N/A' ;
        MOCA = 'N/A' ; 
    end
    
       
    % [ CORRELACIONES ] 
    ROI = strcat( 'ROISignals_' , subID , '.mat' ) ;
    p_ROI = strcat(working_dir, signals_dir , ROI) ;
    S = load (p_ROI) ;
    S = S.ROISignals ; % señales BOLD   
    corr_S = corr(S) ; % cross corr
    hold off ;
    
    % [ SEÑALES BOLD ] 
    for j = 1:116
        sp2 = subplot(2,5,2) ;
        plot (S(:,j)) ;
        title ('Señales BOLD') ;
        xlim( [0 size(S,1)] ) ;
        hold on ;
    end
    hold off ;
    
    
    % [ MATRICES ] 
    sp4 = subplot(2,5,4) ;
    imagesc(corr_S) ;
    title ('matriz correlación') ;
    sp7 = subplot(2,5,7) ;
    spy(corr_S > 0.5) ;
    hold on
    title ('mayor a 0.5') ;
    xlabel(sprintf('', nnz(corr_S > 0.5))) ;
    hold off ;
    
    % [ HISTOGRAMA ] 
    sp3 = subplot(2,5,3) ;
    histogram(corr_S)
    xlim([-1.1 , 1.1]) ;
    title('histograma correlación cruzada') ;
    hold off ;
    
    % [ AUTOCORRELACION ] 
    lags= 75 ;
    for j = 1:116
        y = autocorr ( S(:,j),'NumLags',lags ) ;
        sp5 = subplot(2,5,5) ;
        plot (y) ;
        title ('acf') ;
        hold on ;
    end
    hold off ;
    
    
    % [ MEDIA VS VARIANZA ] 
    m=[] ; v=[] ;
    for j = 1:116
        m(j) = abs(mean ( S(:,j) )) ;
        v(j) = abs(var ( S(:,j) )) ;
    end
    sp6 = subplot(2,5,6) ;
    loglog(m,v,'r*') ;
    xlabel('media') ;
    ylabel('varianza') ;   
    hold off ;
    
    % [ MOVIMIENTO ]
    report_files = dir ( strcat(working_dir,'/RealignParameter/' ,subID)) ;
    for i=3:length(report_files) % busca el archivo
        if contains( string(report_files(i).name) , 'rp')
            nombre = string(report_files(i).name)  ;
            report = load(strcat(working_dir,'/RealignParameter/' ,subID,'/', nombre)) ;
        end
    end
    sp8 = subplot(2,5,8) ;
    plot(report(:,1)) ; 
    hold on ; 
    plot(report(:,2)) ; 
    plot(report(:,3)) ;
    title('movimiento') ;
    legend({'x','y','z'},'Location','southeast','FontSize',6) ;
    xlim([ 0 , size(report,1)] ) ;
    ylabel('mm') ;
    hold off ;
    
    % [ DATOS DEL SUJETO ]
    
    sp9 = subplot(2,5,9) ;
    subID = string(subID) ; fecha = string(fecha) ; imageID = string(imageID) ;
    interlin = 0.15 ; font = 10 ;
    renglon0 = 0.0 ;
    renglon1 = renglon0 + interlin;
    renglon2 = renglon1 + interlin;
    renglon3 = renglon2 + interlin;
    renglon4 = renglon3 + interlin;
    renglon5 = renglon4 + interlin;
    renglon6 = renglon5 + interlin;
    x0 = 0.1 ; x1 = 0.45 ; x2 = 0.7 ; caracteres = 45 ;
    
    text(x0,renglon6, 'Datos del sujeto' ,'FontSize',font+1,'FontWeight','bold') ;
    
    text(x0, renglon5 ,strcat('imageID:',{' '} ,imageID) ,'FontSize',font) ;
    
    text(x0, renglon4 ,estudio,'FontSize',12) ;
    text(x1, renglon4 ,strcat('fecha:',{' '},fecha) ,'FontSize',font) ;
    
    text(x0, renglon3 ,grupo,'FontSize',12) ;
    text(x1, renglon3 ,genero,'FontSize',12) ;
    text(x2, renglon3 ,strcat('edad:',{' '},edad),'FontSize',font) ;
    
    text(x0, renglon2, strcat('MMSE:',{' '}, string(MMSE)) ,'FontSize',font) ;
    text(x1, renglon2, strcat('MOCA:',{' '}, string(MOCA)) ,'FontSize',font) ;

    text(x0, renglon1 ,strcat('Time Points =',TP),'FontSize',font) ;
    text(x2, renglon1, strcat('Quality score:' ,{' '}, string(imgQ)) ,'FontSize',font) ;
    
    if ~ismissing(study_com)
        if strlength(study_com)>=caracteres
            str = char(study_com) ;
            str1 = str(1:caracteres) ;
            corte = find(str1==' ') ;
            corte = corte( size(corte,2) ) ;
            str1 = str(1:corte) ;
            str2 = str(corte+1 : length(str) ) ;
            text(x0, renglon0, str1, 'FontSize',font) ;
            text(x0, renglon0-interlin, str2, 'FontSize',font) ;            
        else
            text(x0, renglon0, study_com, 'FontSize',font) ;            
        end
    end
    axis off ;
    hold off ;
    
    % [ POWER SPECTRUM ]
    %S = S.ROISignals ; % señales BOLD   
    allsignals = [] ;
    signal = [] ;
    TimePoints = size(S,1) ;
    for j=1:116
        signal =  S( (1:TimePoints), j) ;
        signal = normalize(signal);
        signal = reshape ( signal , 1,[] ); % convierto la señal a un vector horizoantal
        allsignals = horzcat( allsignals , signal ) ;
    end
    
    [power , freq ] = pspectrum(allsignals,3) ;
    power = reshape ( power , 1,[] );
    sp10 = subplot(2,5,10) ;
    plot( freq , 10*log10(power)) ;
    title('power spectrum') ;
    xlabel('frecuencias') ;
    hold off ;
    
    % [ TITULO ]
    ss = char(subID) ;
    suptext = strcat( ss(1:3) ,'_', ss(5) ,'_', ss(7:10)) ;
    tlt = suptitle(suptext) ;
    tlt.FontSize = 34 ;
    tlt.Position = [0.165 -0.075 0] ;
    set(tlt,'Interpreter','none') ;
    
    % [ POSICIONES ] 
    % comentar esta parte si no funciona en las proporciones de su pantalla
    % [ left , bottom , width , heigth ]
    % Distance from the left edge 
    % Distance from the bottom edge
    % Distance between the right and left inner edges of the figure
    % Distance between the top and bottom inner edges of the window.
        
    mgleft = 0.035 ; mgcenter = 0.365 ; mgr = 0.74 ; 
    line1 = 0.7 ; line2 = 0.38 ; line3 = 0.05 ;
    w_signals = 0.325 ; h_signals = 0.25 ;
    wfreq = 0.22 ; hfreq = h_signals ;
    wmatrix = 0.25 ; hmatrix = wmatrix ;
    wimg = 0.5 ; himg = wimg ;
    
    pos1 = get(sp1, 'Position')  ; % imagen cerebros 
    new_pos1 = [ -0.085 0.35 wimg himg] ;
    set(sp1, 'Position',new_pos1 ) ; 
    
    pos2 = get(sp2, 'Position') ; % señales BOLD
    new_pos2 = [mgcenter line2 w_signals h_signals] ;
    set(sp2, 'Position',new_pos2 ) ; 
    
    pos3 = get(sp3, 'Position') ; % histograma
    new_pos3 = [mgr line1 wfreq hfreq] ;
    set(sp3, 'Position',new_pos3 ) ; 
    
    pos4 = get(sp4, 'Position') ; % matriz corr color
    new_pos4 = [ mgleft line3 wmatrix/2 hmatrix] ;
    set(sp4, 'Position',new_pos4 ) ; 
    
    pos5 = get(sp5, 'Position') ; % acf
    new_pos5 =[mgcenter line1 w_signals h_signals] ;
    set(sp5, 'Position',new_pos5 ) ;
    
    pos6 = get(sp6, 'Position') ; % media vs varianza
    new_pos6 = [ mgleft+0.19 line2 0.10 0.2] ;
    set(sp6, 'Position',new_pos6 ) ;
    
    pos7 = get(sp7, 'Position') ; % matriz corr spy 
    new_pos7 = [ mgleft+0.092 line3 wmatrix hmatrix] ;
    set(sp7, 'Position',new_pos7 ) ;
    
    pos8 = get(sp8, 'Position') ; % movimiento
    new_pos8 = [mgcenter line3 w_signals h_signals] ;
    set(sp8, 'Position',new_pos8 ) ;
    
    pos9 = get(sp9, 'Position') ; % datos del sujeto
    new_pos9 = [mgr-0.025 line3 wfreq hfreq] ;
    set(sp9, 'Position',new_pos9 ) ;
    
    pos10 = get(sp10, 'Position') ; % power spectrum
    new_pos10 = [mgr line2 wfreq hfreq] ;
    set(sp10, 'Position',new_pos10 ) ;
    
    
    
    % [ GUARDAR ARCHIVOS ] 
    if save==1
        destination = strcat(working_dir,output) ; 
        if ~exist(destination, 'dir')
           mkdir(destination)
        end
        cd(destination) 
        nombre = strcat(subID) ;
        saveas(fig,nombre,'epsc')
    end
    
    
    
   
end 
