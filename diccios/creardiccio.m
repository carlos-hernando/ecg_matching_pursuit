function creardiccio(tipo,delta,limconvini,limconvfin,izmax,demax,segmento)

% Funcion que genera los atomos de un diccionario segun los diferentes
% valores de entrada. Está basada en funciones de Lisandro Lovisolo 
% (triangulos) y Manuel Blanco (convoluciones). Solo se generan los atomos 
% con su valor inicial en la posicion 1 (no a lo largo de todo el segmento
% como con la funcion original "creardiccio.m"), para luego aplicar el MP 
% mediante convolucones circulares de estos atomos.
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
%    Parámetros de salida: Ficheros con estas variables:
%       atomos: matriz con los atomos del diccionario como columnas, por 
%       tanto tiene tamaño -> segmento x #atomos
%       flipatomos: los atomos invertidos en modo espejo para usarlos con
%       las convoluciones circulares del MP
%       parametros: matriz con los datos que definen el atomo, si es
%       triangular: distancia izquierda y distancia derecha; si es
%       convolucion: indice izquierdo e indice derecho
%
%   NOTA: en comparación con el original "creardiccio.m" los parámetros 
%   "centro" desaparecen ya que todos empiezan en la primera posición 
%   del segmento
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 2 abril 2015           |
%             | AUTOR: Carlos Hernando                       |
%             |______________________________________________|
%

diractual = pwd;

%Limite de atomos por matriz para que no haya problemas de memoria en
%Matlab:
disp('Limite de atomos por archivo: ');
limite=2^(20-log2(segmento))
%25
atomos = zeros(segmento,limite);
flipatomos = zeros(segmento, limite);

indice=1;

%altura triangulos:
h=1;

m=0;
N=segmento;

%Creamos una carpeta con el nombre del diccionario: dic+parametros de
%entrada y nos situamos dentro de la carpeta
if tipo==0
    feval('mkdir',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(izmax) ...
    '_' num2str(demax) '_' num2str(segmento)]);
    feval('cd',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(izmax) ...
    '_' num2str(demax) '_' num2str(segmento)]);
else
     feval('mkdir',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(limconvini) '_' num2str(limconvfin) '_' num2str(izmax) ...
    '_' num2str(demax) '_' num2str(segmento)]);
    feval('cd',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(limconvini) '_' num2str(limconvfin) '_' num2str(izmax) ...
    '_' num2str(demax) '_' num2str(segmento)]);
end

diccio = struct('atomos',{},'param',{});

if tipo==0
        
    if delta==0
        
        %Triangular, sin delta     
            
        for j=0:izmax
                        
            for k=0:demax
                
                % Condiciones para evitar deltas, que tambien 
                % se consiguen con longitud de lado 1 
                        
                if (j>0 || k>0) && j~=1 && k~=1 && (k+j+1<=N)
                                        
                    %Asumimos altura=1
                    preatomo=creatriang(j,k,1);
                    
                    m=m+1;
                    atomos(:,m)=zeros(N,1);
                    flipatomos(:,m)=zeros(N,1);
                    atomos(1:1+j+k,m)=preatomo';
                    %En la practica el lado de longitud 0 no existe
                    if j==0 && k==0
                        parametros(:,m)=[1 1]';
                    elseif j==0
                        parametros(:,m)=[1 k]';
                    elseif k==0
                        parametros(:,m)=[j 1]';
                    else
                        parametros(:,m)=[j k]';
                    end
                        
                    %[lado_izq lado_derecho]. El centro es siempre 1+j
                    atomos(:,m)=atomos(:,m)/norm(atomos(:,m));
                    flipatomos(:,m)=circshift(flipud(atomos(:,m)),1);

                    if m==limite
                        diccio(1).atomos=atomos;
                        diccio(1).flipatomos=flipatomos;
                        diccio(1).param=parametros;
                        feval('save',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(izmax) ...
                        '_' num2str(demax) '_' num2str(segmento) '_' num2str(indice) '.mat'], 'diccio');
                        clear diccio;
                        disp('Creando archivo numero: ');
                        indice=indice+1
                        parametros=[];
                        m=0;
                    end
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
                   
                   %Asumimos altura=1
                   preatomo=creatriang(j,k,1);

                   m=m+1;
                   atomos(:,m)=zeros(N,1);
                   flipatomos(:,m)=zeros(N,1);
                   atomos(1:1+j+k,m)=preatomo';
                   %En la practica el lado de longitud 0 no existe
                    if j==0 && k==0
                        parametros(:,m)=[1 1]';
                    elseif j==0
                        parametros(:,m)=[1 k]';
                    elseif k==0
                        parametros(:,m)=[j 1]';
                    else
                        parametros(:,m)=[j k]';
                    end
                   
                   %[lado_izq lado_derecho]. El centro es siempre 1+j
                   atomos(:,m)=atomos(:,m)/norm(atomos(:,m));
                   flipatomos(:,m)=circshift(flipud(atomos(:,m)),1);

                   if m==limite
                        diccio(1).atomos=atomos;
                        diccio(1).flipatomos=flipatomos;
                        diccio(1).param=parametros;
                        feval('save',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(izmax) ...
                        '_' num2str(demax) '_' num2str(segmento) '_' num2str(indice) '.mat'], 'diccio');
                        clear diccio;
                        disp('Creando archivo numero: ');
                        indice=indice+1
                        parametros=[];
                        m=0;
                   end
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
                    
                    [preatomos, param]=Atom_design_new(j,k,limconvini,limconvfin);
                    [A,B]=size(preatomos);
                    
                    for p=1:A
                            
                            m=m+1;
                            atomos(:,m)=zeros(N, 1);
                            flipatomos(:,m)=zeros(N, 1);
                            atomos(1:1+j+k,m)=preatomos(p,:)';
                            % Se guardan los identificadores de cada medio
                            % atomo segun el orden de creación.
                            parametros(1,m) = paramtoindice(param(p,1:2));
                            parametros(2,m) = paramtoindice(param(p,3:4));
                            %Normalizamos los atomos
                            atomos(:,m)=atomos(:,m)/norm(atomos(:,m));
                            flipatomos(:,m)=circshift(flipud(atomos(:,m)),1);                            
                            
                            if m==limite
                                    diccio(1).atomos=atomos;
                                    diccio(1).flipatomos=flipatomos;
                                    diccio(1).param=parametros;
                                    feval('save',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(limconvini) '_' num2str(limconvfin) '_' num2str(izmax) ...
                                    '_' num2str(demax) '_' num2str(segmento) '_' num2str(indice) '.mat'], 'diccio');
                                    clear diccio;
                                    disp('Creando archivo numero: ');
                                    indice=indice+1
                                    parametros=[];
                                    m=0;
                            end
                            
                    end

                end

            end

        end

    else
        
        %Convolucion y con delta
        
        for j=0:izmax
            
            for k=0:demax
                      
               if j~=1 && k~=1 && (k+j+1<=N)
                
                    [preatomos, param]=Atom_design_new(j,k,limconvini,limconvfin);
                    [A,B]=size(preatomos);

                    for p=1:A

                       
                            m=m+1;
                            atomos(:,m)=zeros(N, 1);
                            flipatomos(:,m)=zeros(N, 1);
                            atomos(1:1+j+k,m)=preatomos(p,:)';
                            % Se guardan los identificadores de cada medio
                            % atomo segun el orden de creación.
                            parametros(1,m) = paramtoindice(param(p,1:2));
                            parametros(2,m) = paramtoindice(param(p,3:4));
                            %Normalizamos los atomos
                            atomos(:,m)=atomos(:,m)/norm(atomos(:,m));
                            flipatomos(:,m)=circshift(flipud(atomos(:,m)),1); 
                           
                           if m==limite
                                    diccio(1).atomos=atomos;
                                    diccio(1).flipatomos=flipatomos;
                                    diccio(1).param=parametros;
                                    feval('save',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(limconvini) '_' num2str(limconvfin) '_' num2str(izmax) ...
                                    '_' num2str(demax) '_' num2str(segmento) '_' num2str(indice) '.mat'], 'diccio');
                                    clear diccio;
                                    disp('Creando archivo numero: ');
                                    indice=indice+1
                                    parametros=[];
                                    m=0;
                           end

                    end
                    
               end
               
            end
            
        end
        
    end
     
end

%Guardamos el ultimo si está incompleto
if m~=0
    if tipo==0
        diccio(1).atomos=atomos(:,1:m);
        diccio(1).flipatomos=flipatomos(:,1:m);
        diccio(1).param=parametros;
        feval('save',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(izmax) ...
        '_' num2str(demax) '_' num2str(segmento) '_' num2str(indice) '.mat'], 'diccio');
        clear diccio;
    else
        diccio(1).atomos=atomos(:,1:m);
        diccio(1).flipatomos=flipatomos(:,1:m);
        diccio(1).param=parametros;
        feval('save',['dic_' num2str(tipo) '_' num2str(delta) '_' num2str(limconvini) '_' num2str(limconvfin) '_' num2str(izmax) ...
        '_' num2str(demax) '_' num2str(segmento) '_' num2str(indice) '.mat'], 'diccio');
        clear diccio;
    end
end

cd(diractual);
  