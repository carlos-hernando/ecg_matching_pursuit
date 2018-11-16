function numatomos=tamdiccio(tipo,delta,limconvini,limconvfin,izmax,demax,segmento)

% Funcion que indica el tamaño de un diccionario creado con la 
% funcion "creardiccio_nuevo", es decir, con los diferentes atomos
% situados al principio del segmento, segun los diferentes
% valores de entrada. Está basada en funciones de Lisandro Lovisolo 
% (triangulos) y Manuel Blanco (convoluciones)
% 
%    Parámetros de entrada:
%       tipo: (0)triangulos (1) convoluciones
%       limconvini: numero minimo de convoluciones que se permiten para crear
%       los atomos (solo aplica si tipo=1), si limiteconvini=0 no hay limite inferior.
%       limconvfin: numero maximo de convoluciones que se permiten para crear
%       los atomos (solo aplica si tipo=1), si limiteconvfin=0 no hay
%       limite superior.
%       delta: si se permite la funcion delta (1) o no(0)
%       izmax: distancia maxima entre el comienzo y el centro del atomo
%       demax: distancia maxima entre el centro y el final del atomo
%       segmento: tamaño del segmento en que se divide el ECG
%
%    Parámetros de salida:
%        numatomos: numero de atomos del diccionario
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 19 de octubre 2014     |
%             | AUTOR: Carlos Hernando                       |
%             |______________________________________________|
%

numatomos=0;
N=segmento;

if tipo==0
        
    if delta==0
        
        %Triangular, sin delta
        for j=0:izmax
                                    
            for k=0:demax
                
                % Condiciones para evitar deltas, que tambien 
                % se consiguen con longitud de lado 1 
                        
                if (j>0 || k>0) && j~=1 && k~=1 && (k+j+1<=N)
                    
                    numatomos=numatomos+1;

                end

            end

        end
               
    else

        %Triangular y con delta
        for j=0:izmax
            
            for k=0:demax
                
                % Condiciones para evitar deltas, que tambien 
                % se consiguen con longitud de lado 1 
                        
                if j~=1 && k~=1 && (k+j+1<=N)
                    
                    numatomos=numatomos+1;
                
                end          

            end

        end

    end
    
else
        
    if delta==0
        
        %Convolucion y sin delta
        
        for j=0:izmax
            
            for k=0:demax
                
                if (j>0 || k>0) && j~=1 && k~=1 && (k+j+1<=N)
                    
                    A=Atom_tam_new(j,k,limconvini,limconvfin);
                    
                    for p=1:A
                        
                        numatomos=numatomos+1;

                    end

                end

            end

        end

    else
        
        %Convolucion y con delta
        
        for j=0:izmax
                        
            for k=0:demax
                                      
                if j~=1 && k~=1 && (k+j+1<=N)
                    
                    A=Atom_tam_new(j,k,limconvini,limconvfin);
                    
                    for p=1:A
                        
                        numatomos=numatomos+1;

                    end
                        
                end
                
            end
            
        end
        
    end
        
end