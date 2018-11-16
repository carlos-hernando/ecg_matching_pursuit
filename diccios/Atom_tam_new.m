function numatomos=Atom_tam_new(TL,TR,limconvini,limconvfin)
% Funcion creada a partir de Atom_design_new que indica el numero de 
% atomos que se generarían para unos valores determiandos de TL y TR
% limconvini: limite inferior de # de convoluciones (si ==0 no hay limite 
% inferior)
% limconvini: limite superior de # de convoluciones (si ==0 no hay limite
% superior)
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 26 de noviembre de 2013|
%             | AUTOR: Carlos Hernando Ramiro                |
%             |______________________________________________|
%
%

    %Caso particular de la delta:
    if (TR + TL) == 0
        
        numatomos = 1;
        
    else

        numatomos=0;
    
        [Left_atom, Left_stage]=Half_atom_design_new(TL);
        [Right_atom, Right_stage]=Half_atom_design_new(TR);

        %Para que muestre los parametros en filas y no en columnas
        Left_stage=Left_stage';
        Right_stage=Right_stage';

        if (limconvini > 0) || (limconvfin >0)

            %Hay que guardar solo los atomos que cumplen la condicion
            Left_atompre=Left_atom;
            Right_atompre=Right_atom;
            Left_stagepre=Left_stage;
            Right_stagepre=Right_stage;

            Left_stage=[];
            Right_stage=[];
            Left_atom=[];
            Right_atom=[];

            if limconvini == 0

                contador=1;

                for i=1:size(Left_stagepre,1)
                    if Left_stagepre(i,2)<=limconvfin
                        Left_stage(contador,:)=Left_stagepre(i,:);
                        Left_atom(contador,:)=Left_atompre(i,:);
                        contador=contador+1;
                    end
                end

                contador=1;

                for i=1:size(Right_stagepre,1)
                    if Right_stagepre(i,2)<=limconvfin
                        Right_stage(contador,:)=Right_stagepre(i,:);
                        Right_atom(contador,:)=Right_atompre(i,:);
                        contador=contador+1;
                    end
                end

            elseif limconvfin == 0

                contador=1;

                for i=1:size(Left_stagepre,1)
                    if Left_stagepre(i,2)>=limconvini
                        Left_stage(contador,:)=Left_stagepre(i,:);
                        Left_atom(contador,:)=Left_atompre(i,:);
                        contador=contador+1;
                    end
                end

                contador=1;

                for i=1:size(Right_stagepre,1)
                    if Right_stagepre(i,2)>=limconvini
                        Right_stage(contador,:)=Right_stagepre(i,:);
                        Right_atom(contador,:)=Right_atompre(i,:);
                        contador=contador+1;
                    end
                end

            else

                contador=1;

                for i=1:size(Left_stagepre,1)
                    if ((Left_stagepre(i,2)<=limconvfin) && (Left_stagepre(i,2)>=limconvini))
                        Left_stage(contador,:)=Left_stagepre(i,:);
                        Left_atom(contador,:)=Left_atompre(i,:);
                        contador=contador+1;
                    end
                end

                contador=1;

                for i=1:size(Right_stagepre,1)
                    if ((Right_stagepre(i,2)<=limconvfin) && (Right_stagepre(i,2)>=limconvini))
                        Right_stage(contador,:)=Right_stagepre(i,:);
                        Right_atom(contador,:)=Right_atompre(i,:);
                        contador=contador+1;
                    end
                end
            end
        end 


        lon1=size(Left_stage,1);
        lon2=size(Right_stage,1);
    
        if lon1==0 && TL==0 && lon2 > 0
            for k2=1:lon2
                numatomos=numatomos+1;
            end
        elseif lon2==0 && TR==0 && lon1 > 0
            for k1=1:lon1
                numatomos=numatomos+1;
            end
        elseif lon1 > 0 && lon2 > 0
            for k1=1:lon1
                for k2=1:lon2
                    numatomos=numatomos+1;
                end
            end
        else
            Atoms = [];
            Stages = [];
        end
    
    end
      
        
   
    
   
   
    