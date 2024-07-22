# """export data in other formats apart from julia files"""function [ output_args ] = WriteBrukerRFPulseShapeLau(RF_Shape)


#     TEXT= ['##TITLE= /opt/PV5.1/exp/stan/nmr/lists/wave/',RF_Shape.Name,'\n',...
#     '##JCAMP-DX= 5.00 Bruker JCAMP library\n',...
#     '##DATA TYPE= Shape Data\n',...
#     '##ORIGIN= Technical University of Munich\n',...
#     '##OWNER= <Amanda Nicotina>\n',...
#     '##MISC. ',sprintf('%s',RF_Shape.info),'\n',...
#     '##DATE= ',datestr(now, 'dd/mm/yyyy'),'\n',...
#     '##TIME= ',datestr(now, 'HH:MM:SS'),'\n',...
#     '##MINX= ',num2str(min(RF_Shape.amplitude),'%10.6e'),'\n',...
#     '##MAXX= ',num2str(max(RF_Shape.amplitude),'%10.6e'),'\n',...
#     '##MINY= ',num2str(min(RF_Shape.phase),'%10.5e'),'\n',...
#     '##MAXY= ',num2str(max(RF_Shape.phase),'%10.6e'),'\n',...
#     '##$SHAPE_EXMODE=',RF_Shape.EXMODE,'\n',...
#     '##$SHAPE_TOTROT=',num2str(RF_Shape.TOTROT,'%10.6e'),'\n',...
#     '##$SHAPE_BWFAC=',num2str(RF_Shape.BWFac),'\n',...
#     '##$SHAPE_INTEGFAC=',num2str(RF_Shape.IF),'\n',...
#     '##$SHAPE_REPHFAC=',num2str(RF_Shape.rephasingfactor),'\n',...
#     '##$SHAPE_TYPE=' sprintf('%s',RF_Shape.type),'\n',...
#     '##$SHAPE_MODE= 0\n',...
#     '##MAX_B1_microT=', sprintf('%0.3f',RF_Shape.B1MaxPower),'\n',...
#     '##NPOINTS=',num2str(RF_Shape.length),'\n',...
#     '##XYPOINTS=(XY..XY)','\n'];
    
    
#     FileToWrite=fopen(RF_Shape.Name2save,'w');
#     fprintf(FileToWrite,TEXT);
    
#     for boucle = 1:RF_Shape.length
#         fprintf(FileToWrite,[num2str(RF_Shape.amplitude(boucle),'%10.6e'),', ',num2str(RF_Shape.phase(boucle),'%10.6e'),'\n']);
#     end
#     fprintf(FileToWrite,'##END\n');
#     fclose(FileToWrite);
#     output_args=true;
    
#     % figure
#     % plot(RF_Shape.amplitude)
#     % hold on
#     % plot(RF_Shape.phase,'r')
#     % title('RF_Shape.Name2save')
    
    
#     end
    
    # function exportBrucker(opt,rf,oname)

    #     RF_Shape.info   = '';
    #     RF.type         = 'conventional' ;
        
        
    #     % %%%%%%%%%%%% BandWidth Facor %%%%%%%%%%%%%%%%%%%%%
    #     % ns=opt.N;
    #     % Tp=opt.tf(min(numel(opt.tf),opt.grapeInstance));
    #     % 
    #     % [BW, rffft, f]=bandWidth(ns,rf,Tp,0);
    #     % pause(0.1);
    #     % disp(['BW=',num2str(BW),' Hz'])
        
    #     RF_Shape.BWFac=opt.bandwidthHz(min(opt.grapeInstance,numel(opt.bandwidthHz)))*opt.tf;
        
    #     %%%%%%%%%% Integral Factor %%%%%%%%%%%%%%%%%%%%%%%%%%
    #     % integFac=abs(mean(rf/max(rf)));
    #     % % disp(['integFac=',num2str(integFac)])
    #     % 
    #     % Pout=7.3+20*log10(integFac/.1794*Tp/2e-3);
    #     % disp(['Pout=',num2str(Pout),' dB'])
    #     RF_Shape.IF = trapz(real(rf)/max(real(rf)))/opt.N ;  
        
    #     %%%%%%%%%% Rephasing Factor %%%%%%%%%%%%%%%%%%%%%%%%%
    #     if isfield(opt,'rephaseFactor')
    #         RF_Shape.rephasingfactor = opt.rephaseFactor ;
    #     else
    #         RF_Shape.rephasingfactor=0.0;
    #     end
        
        
    #     %%%%%%%%%%%%%%%%% write pulse %%%%%%%%%%%%%%%%%%%%%%%
    #     RF_Shape.type='ex';
    #     RF_Shape.Filter='ls';
        
    #     if (strcmp(RF_Shape.type,'sat') || strcmp(RF_Shape.type,'st') || strcmp(RF_Shape.type,'ex'))
    #         RF_Shape.TOTROT=90.0;
    #         RF_Shape.EXMODE= 'Excitation';
    #         RF_Shape.NAME_EXT='.exc';
    #     elseif strcmp(RF_Shape.type,'inv')
    #         RF_Shape.TOTROT=180.0;
    #         RF_Shape.EXMODE= 'Inversion';
    #         RF_Shape.NAME_EXT='.inv';
    #         RF_Shape.rephasingfactor=0;
    #     elseif strcmp(RF_Shape.type,'se')
    #         RF_Shape.TOTROT=180.0;
    #         RF_Shape.EXMODE= 'Refocussing';
    #         RF_Shape.NAME_EXT='.rfc';
    #         RF_Shape.rephasingfactor=0;
    #     end
        
        
    #     RF_Shape.insliceripples = 0.01;
    #     RF_Shape.outsliceripples = 0.01;
        
    #     TEMP=rf;
    #     RF_Shape.amplitude=abs(TEMP)/max(abs(TEMP))*100.0;
    #     RF_Shape.phase=angle(TEMP)/pi*180.0-0;
    #     for ii=1:length(RF_Shape.phase)
    #         if RF_Shape.phase(ii)<0
    #             RF_Shape.phase(ii)=360+RF_Shape.phase(ii);
    #         end
    #     end
    #     RF_Shape.length=length(rf);
    #     RF_Shape.integralfactor=abs(mean(TEMP/max(TEMP)));
    #     % RF_Shape.bandwidthfactor=BWFac;
    #     RF_Shape.Name=oname;
    #     RF_Shape.B1MaxPower = max(abs(rf))/267.513 ;
    #     % RF_Shape.Folder='.';
        
        
    #     RF_Shape.Name2save=[oname,'.exc'];
        
    #     % figure
    #     % plot(RF_Shape.amplitude)
    #     % hold on
    #     % plot(RF_Shape.phase,'r')
    #     % plot(RF_Shape.amplitude.*cos(RF_Shape.phase*pi/180),'black')
    #     % plot(RF_Shape.amplitude.*sin(RF_Shape.phase*pi/180),'g')
        
    #     WriteBrukerRFPulseShapeLau(RF_Shape);