function [Atoms, Stages]=Atom_design_new(TL,TR,limconvini,limconvfin)
% Funcion igual que Atom_design de Manuel Blanco pero que devuelve 
% tambien los parametros de convolucion y permite medios atomos TL,TR=0
% El formato de Stages es: [Linileft, #convleft, Liniright, #convright]
% limconvini: limite inferior de # de convoluciones (si ==0 no hay limite 
% inferior)
% limconvini: limite superior de # de convoluciones (si ==0 no hay limite
% superior)
%
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 26 de noviembre de 2013|
%             | AUTOR: Carlos Hernando Ramiro                |
%             |______________________________________________|
%
%

    %Caso particular de la delta:
    if (TR + TL) == 0
        
        Atoms=1;
        Stages=[0 0 0 0];
        
    else %Resto de casos
        
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

        Atoms=zeros(size(Left_stage,1)*size(Right_stage,1),TL+TR+1); %we add +1 to include both limits of the atom the 0 and the 1
        k=1;

        lon1=size(Left_stage,1);
        lon2=size(Right_stage,1);

        
        if lon1==0 && TL==0 && lon2 > 0
            for k2=1:lon2
                Atoms(k2,:)=fliplr(Right_atom(k2,:));
                Stages(k,:)=[0 0 Right_stage(k2,:)];
                k=k+1;
            end
        elseif lon2==0 && TR==0 && lon1 > 0
            for k1=1:lon1
                Atoms(k1,:)=Left_atom(k1,:);
                Stages(k,:)=[Left_stage(k1,:) 0 0];
                k=k+1;
            end
        elseif lon1 > 0 && lon2 > 0
            for k1=1:lon1
                for k2=1:lon2
                    Atoms((k1-1)*size(Right_atom,1)+k2,:)=[Left_atom(k1,:) fliplr( Right_atom(k2,1:TR))];
                    Stages(k,:)=[Left_stage(k1,:) Right_stage(k2,:)];
                    k=k+1;
                end
            end
        else
            Atoms = [];
            Stages = [];
        end
       
    end





