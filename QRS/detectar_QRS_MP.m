function[FNb, FPb, TPb, Seb, Ppb, posQRS, posTP] = detectar_QRS_MP(ficherosenal, diccio, idcal, segini, segfin, itfin, frec, ganancia, senalref, rutabd, tinicomp, sufijo)

% Función que detecta complejos QRS en señales comprimidas mediante el
% método de Matching Pursuit
%
%    Parámetros de entrada:
%       ficherosenal: fichero .mat que guarda los resultados de compresión
%       diccio: carpeta con el diccionario a utilizar
%       idcal: indice de calidad de los resultados de compresion que se
%              quieren emplear de todos los disponibles en el fichero
%       segini: segmento en el que se empieza el analisis
%       segfin: segmento en el que se termina el analisis
%       itfin: iteracion de MP en la que se termina el analisis
%       frec: frecuencia de muestreo del ECG
%       ganancia: ganancia de la conversión ADC de adquisición
%       senalref: señal empleada para caracterizar los latidos candidatos,
%       0: ECG original, 1: ECG recuperada, otro: ECG recuperado solo con
%       los candidatos a latido incluyendo los alternativos
%       rutadb: ruta absoluta a la base de datos a la que pertenece la senal
%       tinicomp: tiempo de inicio de la comprobación de latidos, 0:
%       inicio, otro: tras 5 minutos como dice el estándar
%       sufijo: externsion usada para guardar las detecciones
%
%    Parámetros de salida:
%       FN: numero de falsos negativos
%       FP: numero de falsos positivos
%       TP: numero de latidos correctos
%       Se: sensibilidad de deteccion TP/(TP+FN)
%       Pp: predictibilidad positiva TP/(TP+FP)
%       posQRS: posiciones de todos los latidos detectados
%       posTP: posiciones de todos los aciertos
%       
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 2 abril 2017           |
%             | AUTOR: Carlos Hernando Ramiro                |
%             |______________________________________________|
%
%

%Guardamos el directorio actual para recuperarlo al final
currentfolder = pwd;

%Umbrales configurables
%Umbral para determinar los atomos candidatos a latido
%Valores optimos: segmento 1024 -> 0.4; segmento 2048 -> 0.3
umbcandmax = 0.4;
umbcandmin = umbcandmax / 2;
%Umbral para determinar si un pico es suficientemente alto
umbpic = 0.4;
%Umbral para determinar si una pendiente es suficientemente inclinada
umbpend = 0.4;
%Parac calcular el mínimo de pendiente consideramos como limite un complejo
%QRS de duración tipica 60ms y amplitud maxima 0.15mV, segun estandar EC13
%de la AAMI
minpend = 0.15*ganancia/(0.03*360);
%Umbral para determinar si hay un intervalo sin latido
maxt = 1.5;
%Duracion normal de un complejo onda QRS es de 60ms, entonces cogemos como 
%intervalos para calcular pendientes un cuarto de muestras
intpend = round(0.06 * frec / 4);
%Umbral para determinar si se llega al punto de inflexión
umbinflex = 0.35;
%Divisor desde pico para empezar a calcular pendientes evitando minimos
%locales cercanos al pico
divpend = 3;

%Cargamos los resultados de compresión
feval('load', ficherosenal);

%Tomamos los coeficientes
coef = resultado(idcal).coeficientes;

%Tomamos lo atomos
atomos = resultado(idcal).atomos;

%Obtener todos los ficheros de la carpeta diccio
ficheros=dir(diccio);
cd(diccio);
%Obtener el tamaño del segmento 
%Los dos primeros ficheros son . y ..
feval('load',ficheros(3).name);
tamseg=size(diccio.atomos,1);
clear diccio;

%Obtenemos los indices de los comienzos de coeficientes y atomos
iniparam = [0 find(coef==0)];

%Mínima separacion entre latidos en muestras, 
%asumiendo 300 latidos por minuto (según Pan85)
minsep = round(60 / 300 * frec);

%Maxima separacion entre latidos en muestras,
%asumiendo 30 latidos por minuto
maxsep = round(60 / 30 * frec);

%Duracion estandar de un complejo QRS en muestras (según Guti15)
durqrs = round(60 / 60 * frec);

%Usamos una ventana de verificación de 150 ms de lado, es decir la
%distancia entre anotaciones no puede ser mayor que 150ms: 
%https://www.physionet.org/physiotools/wag/bxb-1.htm
mediavent = floor(0.150 * frec);
%También para calcular la amplitud y la pendiente
busquedapico = round(0.150 * frec);

%Usamos un octavo de la ventana de verificacion (32.5 ms) para buscar 
%el máximo exacto alrededor del pico del atomo
busquedamax = round(0.150 * frec / 8);

%Y la mitad de la ventana para comprobar el ruido
busquedaruido = round(0.150 * frec / 2);

%Tiempos entre latidos de la iteracion anterior y la actual
tentreallant = [];
tentreallact = [];

%Por si se descarta luego el ultimo latido y hay que quitar
%tambien la posible anotacion detectada con él
ultlatdet = false;

%Almacenamos los resultados de todos los segmentos
anotaciondetectada = [];
latidocorrecto = [];
errorlat = [];
FN = 0;
FP = 0;
TP = 0;

%Si ponemos valor 0 de entrada a segfin, cogemos todos los segmentos
if(segfin == 0)
    segfin = length(iniparam) - 1;
end

%Obtenemos las anotaciones de toda la señal
ann = readannot([rutabd '\' resultado(idcal).senal '.atr']);

%Detectamos las anotaciones que corresponden a latidos
%N		Normal beat (displayed as "·" by the PhysioBank ATM and LightWAVE)
%L		Left bundle branch block beat
%R		Right bundle branch block beat
%B		Bundle branch block beat (unspecified)
%A		Atrial premature beat
%a		Aberrated atrial premature beat
%J		Nodal (junctional) premature beat
%S		Supraventricular premature or ectopic beat (atrial or nodal)
%V		Premature ventricular contraction
%r		R-on-T premature ventricular contraction
%F		Fusion of ventricular and normal beat
%e		Atrial escape beat
%j		Nodal (junctional) escape beat
%n		Supraventricular escape beat (atrial or nodal)
%E		Ventricular escape beat
%/		Paced beat
%f		Fusion of paced and normal beat
%Q		Unclassifiable beat
%?		Beat not classified during learning
indlat = find(ann.anntyp=='N' | ann.anntyp=='L' | ann.anntyp=='R' |...
              ann.anntyp=='B' | ann.anntyp=='A' | ann.anntyp=='a' |...
              ann.anntyp=='J' | ann.anntyp=='S' | ann.anntyp=='V' |...
              ann.anntyp=='r' | ann.anntyp=='F' | ann.anntyp=='e' |...
              ann.anntyp=='j' | ann.anntyp=='n' | ann.anntyp=='E' |...
              ann.anntyp=='/' | ann.anntyp=='f' | ann.anntyp=='Q' |...
              ann.anntyp=='?');

%Guardamos los parametros del ultimo latido del segmento anterior
QRSsegant = [];

%Posible ultimo latido del segmento anterior
QRSposant = [];
flagposant = false;
flagincposant = false;
flagyaprocesado = false;

%Señal que se usa para procesar los candidatos
senaldet = [];

%Candiatos que se expanden al segmento siguiente porque no están completos
%en el actual. Distinguimos entre principales y alternativos.
QRScandsigant = [];
QRScandaltsigant = [];

%Si hay que usar el segmento anterior para analizar el candidato
flagsenaldetext = false;

%Si hay que descartar la ultima anotacion en el calculo de los aciertos ya
%que corresponde al segmento siguiente al ultimo y que no se ha procesado
ultimaanotacion = false;

%Posiciones globales de los latdos y los TP
posQRS = [];
posTP = [];

%Analizamos la señal iteracion a iteracion y detectamos los QRS
for seg=segini:1:segfin
    
    %Si en el segmento se ha detectado un latido final parcial
    flagmitadant = false;
    
    %Obtenemos la señal original
    senalorig = resultado(idcal).entrada((seg-1)*tamseg+1:seg*tamseg);
    %Inicializamos la señal a recuperar con todo ceros
    senalrec = zeros(1,tamseg);
    %Tambien generamos la señal recuperada solo con los candidatos
    senalreccand = zeros(1,tamseg);
    %Incluyendo también los alternativos
    senalreccandalt = zeros(1,tamseg);
    %Con los alternativos solo
    senalrecalt = zeros(1,tamseg);
        
    %Inicializamos la lista de QRS
    QRScandini = [];
    k = 0;

    %Iniciamos la lista alternativa de QRS
    QRScandalt = [];
    kalt = 0;
    
    %Candidatos que se expanden al segmento siguiente porque no están completos
    %en el actual. Distinguimos entre principales y alternativos.
    QRScandsig = [];
    QRScandaltsig = [];
    
    %Por deefcto no hay anotacion extra del sgmento siguiente
    anotextra = false;

    %Obtenemos los coeficientes y atomos de un segmento
    coefseg = coef(iniparam(seg)+1:iniparam(seg+1)-1);
    atomseg = atomos(iniparam(seg)+1:iniparam(seg+1)-1,:);

    %Si ponemos valor 0 de entrada a itfin o no hay suficientes iteraciones para
    %llegar a itfin, cogemos todas las iteraciones
    if(itfin == 0 || itfin > length(coefseg))
        itfinseg = length(coefseg);
    else
        itfinseg = itfin;
    end
    
    for it=1:1:itfinseg

        %Tomamos el atomo de la iteracion a partir del diccionario de MP
        if (it==1 || (atomseg(it,1)~=atomseg(it-1,1)))
                    feval('load',ficheros(atomseg(it,1)).name);
        end
        recatomo = coefseg(it) * circshift(diccio.atomos(:,atomseg(it,2)),atomseg(it,3))';
        
        %Reconstruimos la señal
        senalrec = senalrec + recatomo;
        
        %Detectar los QRS
        %El primero siempre lo consideramos candidato
        if(isempty(QRScandini))
            
            %Tomamos la posición del maximo
            maxval = max(abs(recatomo));
            maxvalpri = maxval;
            %Y añadir la posicion como candidata,
            posicion = find(abs(recatomo) == maxval); 
            k = k + 1;
            QRScandini(k,1) = posicion;
            QRScandini(k,2) = sign(coefseg(it));
            QRScandini(k,3) = seg;
            QRScandini(k,4) = it;

            %Reconstruimos la señal solo con los candidatos
            senalreccand = senalreccand + recatomo;

            %Reconstruimos tambien con los alternativos
            senalreccandalt = senalreccandalt + recatomo;

        else
            
            %Comprobar si el maximo del atomo es mayor que un cierto
            %porcentaje del minimo de los maximos
            posiblemaxval = max(abs(recatomo));
            posibleposicion = find(abs(recatomo) == posiblemaxval);
            if(posiblemaxval > umbcandmax * maxval)
                
                %Y añadir la posicion como candidata
                k = k + 1;
                QRScandini(k,1) = posibleposicion;
                QRScandini(k,2) = sign(coefseg(it));
                QRScandini(k,3) = seg;
                QRScandini(k,4) = it;

                %Reconstruimos la señal solo con los candidatos
                senalreccand = senalreccand + recatomo;

                %Reconstruimos tambien con los alternativos
                senalreccandalt = senalreccandalt + recatomo;

            elseif(posiblemaxval > umbcandmin * maxval)

                %Añadir posicion candidata alternativa
                kalt = kalt + 1;
                QRScandalt(kalt,1) = posibleposicion;
                QRScandalt(kalt,2) = sign(coefseg(it));
                QRScandalt(kalt,3) = seg;
                QRScandalt(kalt,4) = it;

                %Reconstruimos con los candidatos alternativos por si acaso
                senalreccandalt = senalreccandalt + recatomo;
                
                %Reconstruimos solo con los alternativos
                senalrecalt = senalrecalt + recatomo;

            end
            
        end
        
    end
    
    %Segmento anterior para el caso en que sea necesario por si el latido está
    %al principio del segmento e incompleto
    senaldetant = senaldet;

    %Usamos la señal original o la recuperada para el analisis de 
    %deteccion de latidos a partir de los candidatos
    if(senalref == 0)
        senaldet = senalorig;
    elseif (senalref == 1)
        %Tenemos que usar la señal en unidades originales, para poder usar
        %el límite de mínima pendiente (nos faltaría el offset, pero los
        %valores absolutos no importan)
        senaldet = senalrec * resultado(idcal).normas(seg);
    else
        senaldet = senalreccandalt * resultado(idcal).normas(seg);
    end
    
    %Guardamos los latidos
    QRSseg = [];
    
    %Bandera de que es la primera busqueda
    flagini = true;
    tamant = 0;
    
    %Lista de posiciones de maximos de picos
    listaposmaxpico = [];
    
    %Grupo incial de posibles candidatos
    QRScand = QRScandini;
    
    %Bandera para decidir si se hacen reintentos
    flagotra = true;
   
    %Si reintento o deteccion de posible ultimo latido
    while(flagotra || flagposant)
        
        %Ordenamos los posibles latidos en el tiempo
        QRScand = sortrows(QRScand, 1);
        
        %Añadimos los candidatos que provienen del segmento anterior
        if(flagini)
            QRScand = [QRScandsigant' QRScand']';
        end

        %Almacenamos las amplitudes de pico, el valor maximo y 
        %las pendientes de los candidatos
        picos = [];
        pendientes = [];
        listaposmaxpico = [];
        
        p=1;
        while(p<=size(QRScand,1))
            
            %Limites de ventana donde buscamos el pico, es la mitad de los
            %limites de verificacion
            minlim = QRScand(p,1) - busquedamax;
            maxlim = QRScand(p,1) + busquedamax;
            
            %Bajamos la bandera de candidato para el segmento siguiente
            %por defecto
            flagcandsig = false;
            
            %Si nos salimos del segmento para atrás, añadimos el segmento
            %anterior
            if(minlim < 1)
                if(~isempty(senaldetant))
                    flagsenaldetext = true;
                    senaldet = [senaldetant senaldet];
                    minlim = minlim + tamseg;
                    maxlim = maxlim + tamseg;
                else
                    %Si es el primer segmento y no hay anterior, se
                    %evalua en este mismo segmento
                    minlim = 1;
                end
            %Si nos salimos del segmento hacia delante, guardamos el
            %candidato para evaluarlo en el siguiente segmento, pero lo
            %evaluamos aquí de todas formas por si acaso hay datos
            %suficientes
            elseif(maxlim > tamseg)
                
                %Para su evaluacion en este mismo segmento
                maxlim = tamseg;
                
                %Indicamos que ya se añadió el candidato para el siguiente
                %segmento
                flagcandsig = true;
                
                %Si no es reintento lo ponemos en el grupo de principales
                if(flagini)
                
                    QRScandsig = [QRScandsig' QRScand(p,:)']';
                    %Actualizamos el tiempo para ponerlo acorde al segmento
                    %siguiente
                    QRScandsig(end,1) = QRScandsig(end,1) - tamseg; 
                
                %Si es reintento y no se añadio ya a ninguno de los dos grupos
                else
                    
                    if(~isempty(QRScandsig) && isempty(find(QRScandsig(:,1) == QRScand(p,1))) &&...
                       ~isempty(QRScandaltsig) && isempty(find(QRScandaltsig(:,1) == QRScand(p,1))))
                   
                        QRScandaltsig = [QRScandaltsig' QRScand(p,:)']';
                        %Actualizamos el tiempo para ponerlo acorde al segmento
                        %siguiente
                        QRScandaltsig(end,1) = QRScandaltsig(end,1) - tamseg;
                    end

                end
                
            end
                            
            %Si el coeficiente era negativo, analizamos la señal
            %negativa porque el pico es hacia abajo
            if(QRScand(p,2) < 0)
                senaldet = -senaldet;
            end
            
            %Obtenemos el maximo y su posicion
            [maxpico, posmaxpico] = max(senaldet(minlim:maxlim));
            posmaxpicotot = posmaxpico + minlim - 1;
            %La posicion del maximo es siempre relativa al segmento actual
            if(flagsenaldetext)
                posmaxpicotot = posmaxpicotot - tamseg;
            end
            
            %Comprobamos si ese pico ya está como candidato
            picorepe = find(listaposmaxpico == posmaxpicotot);
            
            %Si el pico no está repetido
            if(isempty(picorepe))
                
                %Comprobamos si se corresponde con el posible ultimo latido
                %del segmento anterior, para indicar que ya se procesa aquí
                %y no hay que añadirlo luego
                if(flagincposant && posmaxpicotot==QRSposant(1))
                    flagyaprocesado = true;
                end
                
                %Lo añadimos a la lista de posiciones de picos
                listaposmaxpico=[listaposmaxpico posmaxpicotot];
                %Lo usamos como nueva posicion de candidato
                QRScand(p,1) = posmaxpicotot;
                %Actualizamos los limtes de busqueda
                minlim = QRScand(p,1) - busquedapico;
                maxlim = QRScand(p,1) + busquedapico;
                %Si nos salimos del segmento para atrás, añadimos el segmento
                %anterior
                if(minlim < 1)
                    
                    if(~isempty(senaldetant))
                    
                        %Si no habíamos extendido ya el segmento
                        if(flagsenaldetext == false)
                            flagsenaldetext = true;
                            if(QRScand(p,2) < 0)
                                senaldet = [-senaldetant senaldet];
                            else
                                senaldet = [senaldetant senaldet];
                            end
                        end

                        minlim = minlim + tamseg;
                        maxlim = maxlim + tamseg;
                        %Recuperamos la posicio extendida del pico
                        posmaxpicotot = posmaxpicotot + tamseg;
                        
                    else
                        %Si es el primer segmento y no hay anterior, se
                        %evalua en este mismo segmento
                        minlim = 1;
                    end
                    
                elseif(maxlim > tamseg)
                    
                    %Para su evaluacion en este mismo segmento
                    maxlim = tamseg;

                    %Si no se añadio ya antes
                    if(~flagcandsig)
                        
                        %Si no es reintento lo ponemos en el grupo de principales
                        if(flagini)

                            QRScandsig = [QRScandsig' QRScand(p,:)']';
                            %Actualizamos el tiempo para ponerlo acorde al segmento
                            %siguiente
                            QRScandsig(end,1) = QRScandsig(end,1) - tamseg; 

                        %Si es reintento y no se añadio ya a ningun grupo
                        else
                            
                            %Si no está ya entre los principales
                            if(isempty(QRScandsig) ||...
                               (~isempty(QRScandsig) && isempty(find(QRScandsig(:,1) == QRScand(p,1)))))
                              
                                %Y tampoco entre los alternativos
                                if(isempty(QRScandaltsig) ||...
                                 ~isempty(QRScandaltsig) && isempty(find(QRScandaltsig(:,1) == QRScand(p,1))))

                                QRScandaltsig = [QRScandaltsig' QRScand(p,:)']';
                                %Actualizamos el tiempo para ponerlo acorde al segmento
                                %siguiente
                                QRScandaltsig(end,1) = QRScandaltsig(end,1) - tamseg;

                                end
                                
                            end

                        end
                        
                    end 
                    
                end
                   
                %Desde el maximo buscamos los puntos de inflexión
                %Hacia detrás
                    
                nonstop = true;
                q = posmaxpicotot;
                inflexuno = minlim;

                %Obtenemos el valor minimo hacia delante
                minpico = min(senaldet(minlim:q));
                %Distancia entre el valor maximo y el minimo
                distmaxmin = maxpico - minpico;

                %Saltamos del pico hacia abajo para evitar minimos locales
                %cercanos a él como dobles picos
                valmitad = maxpico - (distmaxmin/divpend);
                posmitad = find(senaldet(minlim:q) <= valmitad);
                %Si hay valores, nos colocamos ahí, si no medimos la
                %pendiente desde el principio prouqe no será muy
                %pronunciada
                if(~isempty(posmitad))
                    q = posmitad(end) + minlim - 1;
                end

                %Nos vamos hacia detras para tener muestras delante para
                %calcular la pendiente
                q = q - intpend;

                %Si no hemos llegado al limite seguimos calculando
                %pendientes para buscar el punto de inflexion
                if(q > minlim)

                    %Calculamos la pendiente inicial
                    pendini = (senaldet(q + intpend) - senaldet(q)) / intpend;
                    q = q - 1;

                    %Hasta que lleguemos al limite o la pendiente sea menor que
                    %inicial segun el umbral establecido
                    while(nonstop && q > minlim)

                        %Calculamos la pendiente actual
                        pendactual = (senaldet(q + intpend) - senaldet(q)) / intpend;

                        %Comprobamos si es mucho menor que la inicial
                        if(pendactual < umbinflex * pendini)
                            inflexuno = q + ceil(intpend/2);
                            nonstop = false;
                        end

                        q = q - 1;

                    end

                end
                
                %Hacia delante                     
                nonstop = true;
                q = posmaxpicotot;
                inflexdos = maxlim;

                %Obtenemos el valor minimo hacia delante
                minpico = min(senaldet(q:maxlim));
                %Distancia entre el valor maximo y el minimo
                distmaxmin = maxpico - minpico;

                %Saltamos del pico hacia abajo para evitar minimos locales
                %cercanos a él como dobles picos
                valmitad = maxpico - (distmaxmin/divpend);
                posmitad = find(senaldet(q:maxlim) <= valmitad);
                %Si hay valores, no colocamos ahí, si no, medimos la
                %pendiente desde el principio prouqe no será muy
                %pronunciada
                if(~isempty(posmitad))
                    q = q + posmitad(1) - 1;
                end

                %Nos vamos hacia delante para tener muestras atras para
                %calcular la pendiente
                q = q + intpend;

                %Si no hemos llegado al limite seguimos calculando
                %pendientes para buscar el punto de inflexion
                if(q < maxlim)

                    %Calculamos la pendiente inicial
                    pendini = (maxpico - senaldet(q)) / intpend;
                    q = q + 1;

                    %Hasta que lleguemos al limite o la pendiente sea menor que
                    %inicial segun el umbral establecido
                    while(nonstop && q < maxlim)

                        %Calculamos la pendiente actual
                        pendactual = (senaldet(q - intpend) - senaldet(q))/ intpend;

                        %Comprobamos si es mucho menor que la inicial
                        if(pendactual < umbinflex * pendini)
                            %Lo despalzamos al medio de la recta que
                            %calcula la pendiente
                            inflexdos = q - ceil(intpend/2);
                            nonstop = false;
                        end

                        q = q + 1;

                    end

                end
                
                %Comprobamos que no es ruido. Es ruido si hay picos
                %similares alrededor.
                minlim=max(1,inflexdos-busquedamax);
                maxlim=min(length(senaldet),inflexuno+busquedamax);
                valinfdel = min(senaldet(minlim:inflexdos));
                valinfdet = min(senaldet(inflexuno:maxlim));
                umbraldel = maxpico - ((maxpico - valinfdel)/2);
                umbraldet = maxpico - ((maxpico - valinfdet)/2);
                %Hacia delante
                maxlim = min(length(senaldet), posmaxpicotot+busquedaruido);
                ruidodel = find(senaldet(inflexdos:maxlim) > umbraldel);
                %Hacia detras
                minlim = max(1, posmaxpicotot-busquedaruido);
                ruidodet = find(senaldet(minlim:inflexuno) > umbraldet);
                
                if(isempty(ruidodel) || isempty(ruidodet))
                
                    %Calculamos la diferencia entre el pico y el punto de inflexion más 
                    %alto
                    liminf = max(senaldet(inflexuno),senaldet(inflexdos));
                    picos(p) = maxpico - liminf;

                    %Calculamos ambas pendientes
                    penduno = (maxpico - senaldet(inflexuno)) / (posmaxpicotot - inflexuno);
                    penddos = (maxpico - senaldet(inflexdos)) / (inflexdos - posmaxpicotot);
                    %Y nos quedamos con la más pequeña.
                    pendientes(p) = min(penduno, penddos);
                    
                    %Volvemos a poner la señal original
                    %Si la hemos extendido
                    if(flagsenaldetext)
                        senaldet = senaldet(tamseg+1:end);
                        flagsenaldetext = false;
                    end
                    %Si era negativa
                    if(QRScand(p,2) < 0)
                        senaldet = -senaldet;
                    end
                    
                    %Pasamos al latido siguiente
                    p = p + 1;
                    
                else
                    
                    %Volvemos a poner la señal original
                    %Si la hemos extendido
                    if(flagsenaldetext)
                        senaldet = senaldet(tamseg+1:end);
                        flagsenaldetext = false;
                    end
                    %Si era negativa
                    if(QRScand(p,2) < 0)
                        senaldet = -senaldet;
                    end
                    
                    %Eliminamos el latido actual de la lista de candidatos
                    QRScand(p,:) = [];
                    
                end
                 
            else
               
                %Volvemos a poner la señal original
                %Si la hemos extendido
                if(flagsenaldetext)
                    senaldet = senaldet(tamseg+1:end);
                    flagsenaldetext = false;
                end
                %Si era negativa
                if(QRScand(p,2) < 0)
                    senaldet = -senaldet;
                end
                
                %Eliminamos el latido actual de la lista de candidatos
                QRScand(p,:) = [];
                
            end

        end
        
        %Añadimos el posible ultimo latido del segmento anterior, si
        %estamos reintentando el primer intervalo porque hay demasiado
        %tiempo sin latido y si no lo hemos procesado ya como uno 
        %de los candiatos incompletos que se calculan de nuevo en el
        %segmento siguiente
        if(flagincposant)

            if(~flagyaprocesado)
                            
                QRSnuevo=QRSposant(1:4);
                QRScand = [QRSnuevo' QRScand']';
                picos = [QRSposant(5) picos];
                pendientes = [QRSposant(6) pendientes];
                
            else
                
                flagyaprocesado = false;
            
            end
            
            flagincposant = false;
             
        end
        
        %Añadimos el ultimo latido del segmento anterior a la lista de
        %latidos si es la primera busqueda y actualizamos la lista de
        %candidatos iniciales
        if(flagini)
            QRSseg = [QRSsegant' QRSseg']';
            QRScandini = QRScand;
        end
        
        %Si tenemos candidatos para procesar
        if(~isempty(QRScand))
            
            if(isempty(QRSseg))

                %Tomamos el pico mayor como bueno
                [picomayor, pospicomayor] = max(picos);
                QRSseg(1,:) = QRScand(pospicomayor,:);
                %Añadimos el valor de pico, el de pendiente y el maximo
                QRSseg(1,5) = picos(pospicomayor);
                QRSseg(1,6) = pendientes(pospicomayor);
                %Indice para colocarnos en el latido de referencia
                k = 1;

            else %Si ya hay latidos de antes
                
                %Caso especial en el que no hay latido anterior que usar
                %como referencia y se emplea el posterior
                flagrefpost = false;

                %Tomamos como referencia el latido que marca el limite
                %inicial del intervalo
                %Hay que poner el indice correspondiente dependiendo de si
                %es reintento o no
                if(flagini)
                    ind = 1;
                else
                    ind = find(QRSseg(:,1) == treintento);
                end
                
                %Si no hay un latido previo dsiponible, es decir si es el
                %inicio de la señal, miramos haber si hay al menos
                %uno posterior para usar como referencia,que sería
                %el primero de la lista
                if(isempty(ind) && size(QRSseg,1) > 0)
                    ind = 1;
                    flagrefpost = true;
                end                  
                  
                %Si al final tenemos un indice adecuado
                if(~isempty(ind))
                    
                    %Añadimos sus valores en la primera posicion de picos y pendientes
                    picos = [QRSseg(ind,5) picos];
                    pendientes = [QRSseg(ind,6) pendientes];
                    %Ponemos lo valores correspondientes de pico maximo
                    picomayor = picos(1);
                    pospicomayor = 1;
                    %Lo añadimos tambien al principio de la lista de candidatos
                    %para que los indices con picos y pendientes coincidan
                    QRScand = [QRSseg(ind,1:4)' QRScand']';
                    %Indice para colocarnos en el latido de referencia
                    k = 1;
                   
                %Si no,no se lleva a cabo la busqueda porque no hay
                %referencia
                else
                    
                    picomayor = 0;
                    
                end

            end
            
        else
            
            %Se indica que no hay candidatos disponibles
            picomayor = 0;
            
        end      
        
        if(picomayor > 0)
            
            %Si no es reintento
            if(flagini)
                
                %Evaluamos los posibles picos posteriores al mayor 
                %en comparacion con el anterior
                anterior = pospicomayor;
                for p=pospicomayor+1:1:length(picos)

                    if(pendientes(p) > minpend && picos(p) > (umbpic*picos(anterior)) && pendientes(p) > (umbpend*pendientes(anterior)))

                        %Si hay separacion suficiente
                        if(QRScand(p,1) - QRScand(anterior,1) >= minsep) 

                            %Tomamos los datos del candidato
                            QRSnuevo = QRScand(p,:);
                            %Añadimos el valor de pico y el de pendiente
                            QRSnuevo(5) = picos(p);
                            QRSnuevo(6) = pendientes(p);
                            %Lo añadimos como nuevo latido
                            k = size(QRSseg,1) + 1;
                            QRSseg(k,:) = QRSnuevo;
                            %Y renombramos al anterior
                            anterior = p;
                            
                        %Si tiene mayor pendiente, se sustituye
                        elseif(pendientes(p) > pendientes(anterior))
                            
                            %Tomamos los datos del candidato
                            QRSnuevo = QRScand(p,:);
                            %Añadimos el valor de pico y el de pendiente
                            QRSnuevo(5) = picos(p);
                            QRSnuevo(6) = pendientes(p);
                            %Sustituimos al latido anterior
                            k = size(QRSseg,1);
                            QRSseg(k,:) = QRSnuevo;
                            %Y renombramos al anterior
                            anterior = p;
                            
                        end
                        
                    end

                end
                
                %Evaluamos los posibles picos anteriores al mayor 
                %en comparacion con el posterior
                posterior = pospicomayor;
                for p=pospicomayor-1:-1:1
                   
                    if(pendientes(p) > minpend && picos(p) > (umbpic*picos(posterior)) && pendientes(p) > (umbpend*pendientes(posterior)))

                        
                        %Si tenemos separacion suficente con el posterior,
                        %que es el primero, lo añadimos como nuevo latido
                        if(QRScand(posterior,1)-QRSseg(1,1) >= minsep)
                            
                            %Tomamos los datos del candidato
                            QRSnuevo = QRScand(p,:);
                            %Añadimos el valor de pico y el de pendiente
                            QRSnuevo(5) = picos(p);
                            QRSnuevo(6) = pendientes(p);
                            %Lo añadimos como nuevo latido
                            k = size(QRSseg,1) + 1;
                            QRSseg(k,:) = QRSnuevo;
                            %Y reordenamos
                            QRSseg = sortrows(QRSseg, 1);
                            %Y renombramos al posterior
                            posterior = p;

                         %Si tiene mayor pendiente, se sustituye
                         elseif(pendientes(p) > pendientes(posterior))

                            %Tomamos los datos del candidato
                            QRSnuevo = QRScand(p,:);
                            %Añadimos el valor de pico y el de pendiente
                            QRSnuevo(5) = picos(p);
                            QRSnuevo(6) = pendientes(p);
                            %Sustituimos al latido posterior
                            QRSseg(1,:) = QRSnuevo;
                            %Y renombramos al posterior
                            posterior = p;

                        end
                     
                    end
               
                end
                
            %Si es reintento, evaluamos solo hacia delante y teniendo en
            %cuenta todos los latidos que ya hay detectados en el segmento
            %Comparandonos siempre con el primero y añadiendo solo uno
            %Ademas reducimos los umbrales a la mitad y damos preferencia
            %al que tenga la mayor pendiente
            else  
                
                QRSnuevo=[];
                
                for p=pospicomayor+1:1:length(picos)

                    if(pendientes(p) > minpend && picos(p) > (umbpic*picos(1)/2) && pendientes(p) > (umbpend*pendientes(1)/2))

                       %Si hay separacion suficiente con el primero
                       if(QRScand(p,1) - QRScand(1,1) >= minsep)
                           %Y con el limite posterior si lo hubiere
                           if(flagrefpost || ind + 1 > size(QRSseg,1) ||...                               
                              QRSseg(ind+1,1) - QRScand(p,1) >= minsep)

                               %Si no habia ninguno valido
                               if(isempty(QRSnuevo))
                                   
                                   %Tomamos los datos del candidato
                                   QRSnuevo = QRScand(p,:);
                                   %Añadimos el valor de pico y el de pendiente
                                   QRSnuevo(5) = picos(p);
                                   QRSnuevo(6) = pendientes(p);
                                
                               %Si ya habia algun nuevo valido hay que
                               %comprobar que éste es mejor en terminos de
                               %pendiente
                               elseif(pendientes(p) > QRSnuevo(6))
                                  
                                   %Cambiamos el nuevo latido
                                   QRSnuevo = QRScand(p,:);
                                   %Añadimos el valor de pico y el de pendiente
                                   QRSnuevo(5) = picos(p);
                                   QRSnuevo(6) = pendientes(p);
                                    
                               end
                                
                           end
                           
                       end
                       
                    end
                    
                end
                
                %Si hemos conseguido algún latido nuevo
                if(~isempty(QRSnuevo))
                    %Lo añadimos como posible anterior
                    if(flagposant == true)

                        QRSposant = QRSnuevo;
                        QRSposant(1) = QRSposant(1) - tamseg;

                    else %O como nuevo latido

                        k = size(QRSseg,1) + 1;
                        QRSseg(k,:) = QRSnuevo;
                        %Y reordenamos
                        QRSseg = sortrows(QRSseg, 1);
                        %Cogemos el indice del nuevo pico para colocarnos en él
                        k = find(QRSseg(:,1) == QRScand(pospicomayor,1));

                    end
                    
                end
                
            end
            
        end
            
        if(flagposant == false)
                        
            %Comprobamos los tiempos entre latidos
            tentre = [];
            for n=2:size(QRSseg,1) + 1

                %Guardamos la duracion del intervalo, el tiempo inicial y el
                %final
                if(n == size(QRSseg,1)+1)
                    tentre(n-1,1) = tamseg - 1 - QRSseg(n-1,1);
                    tentre(n-1,2) = QRSseg(n-1,1) + 1;
                    tentre(n-1,3) = tamseg;
                else
                    tentre(n-1,1) = QRSseg(n,1) - QRSseg(n-1,1);
                    tentre(n-1,2) = QRSseg(n-1,1) + 1;
                    tentre(n-1,3) = QRSseg(n,1) - 1;
                end

            end

            %Si no tenemos latido del segmento anterior, nos faltaría el primer
            %intervalo, por lo que lo añadimos
            if(~isempty(QRSseg(1,1)) && QRSseg(1,1) > 0)

                tentreaux = [];
                tentreaux(1) = QRSseg(1,1);
                tentreaux(2) = 1;
                tentreaux(3) = QRSseg(1,1) - 1;
                tentre = [tentreaux' tentre']';

            end

            tamactual = size(tentre,1);

            %Si es la primera vez
            if(flagini)

                %Desactivamos el flagini
                flagini = false;

                %Inicializamos los reintentos a 2, que significa que no se
                %ha reintentado ni con el grupo principal ni con el de
                %candidatos
                reintentar = 2 .* ones(1,tamactual);

            %Si hay intervalos nuevos estaran en el que ya se busco y/o 
            %el nuevo, ambos al final
            elseif(tamant < tamactual)

                noreintentar = find(reintentar == 0);
                
                %El que se reintento y se sacó un nuevo latido hay que
                %comprobarlo de nuevo prou peude seguir siendo
                %demasiado largo pese al nuevo latido
                reintentar = [zeros(1,length(noreintentar)-1) (2 .* ones(1,tamactual - length(noreintentar)+1))];

            end

            %Guardamos el tamaño de tentre para comparar en el futuro
            tamant = size(tentre, 1);

            %Actualizamos la lista de intervalos para esta busqueda
            tentreallact = [tentreallant tentre(:,1)'];

            %Ver cual es el primer intervalo que hay que reintentar
            %para empezar a evaluar desde ahi
            noreintentados = find(reintentar > 0);
            if(isempty(noreintentados))
                n = size(tentre,1) + 1;
            else
                n = noreintentados(1);
            end
            
            %Comparamos si algún intervalo es mucho mayor que el de referencia,
            %porque nos faltarían latidos
            treintento = [];
            flagotra = false;
            while(n <= size(tentre,1) && ~flagotra)

                %El tiempo de referencia es la mediana de los tiempos anteriores,
                %hasta un maximo de 30, excepto si tenemos menos de 4 antes 
                %(el primero no cuenta) que usamos los que tengamos 
                %disponibles menos el mismo
                naux = n;
                if(naux == size(tentre,1))
                    naux = naux - 1;
                end
                if((length(tentreallant) - 1 + naux) < 5)

                    tentreaux = [];
                    %Comprobamos si tenemos hasta 4 intervalos antes,
                    %(el primero y el ultimo no cuentan)
                    %si no cogemos lo maximo que tengamos
                    limitaux = length(tentreallact) - 1;
                    cont = 0;
                    for j=2:limitaux
                        if(cont < 4)
                            %Usamos todos los disponibles menos él mismo
                            if(tentreallact(j) ~= tentre(n,1))
                                tentreaux = [tentreaux tentreallact(j)];
                                cont = cont + 1;
                            end
                        else
                            break;
                        end
                    end

                    if(isempty(tentreaux))
                        tref = tentre(n,1);
                    else
                        tref = median(tentreaux);
                    end

                else

                    interult = length(tentreallant)+ n - 1;

                    %Tenemos en cuenta los 30 ultimos latidos
                    trefliminf = max(2,interult-29);
                    tref = median(tentreallact(trefliminf:interult));

                end

                if(tentre(n,1) > min(tref * maxt, maxsep))

                    %Usamos como nuevo grupo de candidatos el que está en el
                    %intervalo
                    newcand = find(QRScandini(:,1) >= tentre(n,2) & QRScandini(:,1) <= tentre(n,3));
                    QRScand = QRScandini(newcand,:);
                    
                    %Añadimos los incompletos principales del segmento anterior
                    if(reintentar==1)
                        QRScand = [QRScandsigant' QRScand']';
                    end
                    
                    %Si no hay principales o ya se reintentaron
                    if(isempty(QRScand) || reintentar(n) == 1)

                        %Se da por reintentado por completo
                        reintentar(n) = 0;

                        if(~isempty(QRScandalt))

                            newcandalt = find(QRScandalt(:,1) >= tentre(n,2) & QRScandalt(:,1) <= tentre(n,3));
                            QRScand = [QRScand' QRScandalt(newcandalt,:)']';

                        end
                        
                        %Añadimos los incompletos alternativos del segmento anterior
                        if(reintentar==1)
                            QRScand = [QRScandaltsigant' QRScand']';
            
                        end

                    else

                        %Se da por reintentado con el grupo principal
                        reintentar(n) = 1;

                    end

                    %Si hay o si es el primer intervalo y hay un posible ultimo
                    %latido del segmento anterior
                    if(~isempty(QRScand) || (n==1 && ~isempty(QRSposant)))

                        %Indicamos el tiempo de inicio del intervalo de
                        %reintento donde esta el latido referencia
                        treintento = tentre(n,2) - 1;
                        %Y activamos la bandera
                        flagotra = true;

                        %Si es el segundo caso tenemos que activar la
                        %bandera
                        if(n==1 && ~isempty(QRSposant))

                            flagincposant = true;

                        end

                    end

                %Si el intervalo es suficientemente pequeño
                else

                    %Se da por reintentado por completo
                    reintentar(n) = 0;

                end

                n= n + 1;

            end

            %Si ya se han reintentado todos los intervalos, queda conseguir el posible latido
            %del ultimo segmento, porque no se sabe como de largo es el
            %ultimo intervalo y puede que nos lo estemos dejando sin
            %detectar (si no es el ultimo segmento)
            if(~flagotra && ~flagposant && seg<segfin)
                
                %Limpiamos el latido guardado anteriormente
                QRSposant = [];
                %Usamos como nuevo grupo de candidatos el que está en el
                %ultimo intervalo
                n = size(tentre,1);
                newcand = find(QRScandini(:,1) >= tentre(n,2) & QRScandini(:,1) <= tentre(n,3));
                QRScand = QRScandini(newcand,:);
                %Y los alternativos si hay
                if(~isempty(QRScandalt))

                    newcandalt = find(QRScandalt(:,1) >= tentre(n,2) & QRScandalt(:,1) <= tentre(n,3));
                    QRScand = [QRScand' QRScandalt(newcandalt,:)']';

                end

                %Si hay candidatos, se busca el posible ultimo latido
                if(~isempty(QRScand))

                    %Indicamos el tiempo de inicio del intervalo de
                    %reintento
                    treintento = tentre(n,2) - 1;
                    %Y activamos la bandera
                    flagposant = true;

                end

            end
            
        else
            
            %Deshabilitar la bandera porque ya se busco
            flagposant = false;
            
        end
        
    end
    
    %Indicamos si el ultimo latido del segmento anterior se sustituyó
    sustultlat = false;
    %Quitamos el primer latido si es el del segmento anterior
    if(~isempty(QRSsegant))

        QRSseg = QRSseg(2:size(QRSseg,1),:); 

    else
        
        %Marcamos que el ultimo latido fue sustituido
        if(seg > segini)
            sustultlat = true;
        end
        
    end
    
    %Si hay, guardamos el ultimo latido del segmento anterior
    if(size(QRSseg,1) > 0)
        
        QRSsegant = QRSseg(size(QRSseg,1),:);
        %El tiempo del latido lo tenemos que poner correspondiente
        %al segmento siguiente
        QRSsegant(1) = QRSsegant(1) - tamseg;
        
    end
    
    %Guardamos todos los intervalos
    tentreallant = [tentreallant tentre(1:size(tentre,1)-1,1)'];
    
    %Añadimos el offset temporal del segmento comprimido
    %para compararlos con las anotaciones
    tini = (seg-1) * tamseg;
    QRSseg(:,1) = QRSseg(:,1) + tini;
    
    %Obtenemos las anotaciones para ese segmento
    %El primer tiempo que es el de la primera muestra es el indice cero
    tfin = (seg-1) * tamseg + tamseg - 1;
    indicesanotseg = find(ann.time(indlat) >= tini & ann.time(indlat) <= tfin).';
    if(isempty(indicesanotseg))
        indicesanotseg = [];
    end
    
    %Si hay ya alguna anotacion empleada,
    %se tiene en cuenta la ultima del segmento anterior en la comprobacion 
    if(~isempty(anotaciondetectada))
        
        if(~isempty(indicesanotseg))
            
            indicesanotseg = [indicesanotseg(1)-1 indicesanotseg];
            indiceanotsegant = indicesanotseg(length(indicesanotseg));
            
        else
            
            indicesanotseg = indiceanotsegant;
            
        end
        
    end
    
    %Si no es el segmento con la ultima anotacion, añadimos una anotacion
    %mas al final para la comprobacion
    if(indicesanotseg(length(indicesanotseg)) < length(indlat))

        indicesanotseg = [indicesanotseg indicesanotseg(length(indicesanotseg))+1];
        anotextra = true;

        if(seg == segfin)
            ultimaanotacion = true;
        end

    end
    
    %Sumamos uno a las posiciones para convertirlas en indices de Matlab
    lattime = ann.time(indlat(indicesanotseg))' + 1;
    lattipo = ann.anntyp(indlat(indicesanotseg))';
    
    %Comprobar si las anotaciones se detectan y si los latidos son correctos
    inianotseg = indicesanotseg(1) - 1;
    nuevasanot = indicesanotseg(length(indicesanotseg)) - size(anotaciondetectada,2);
    anotaciondetectadaseg = zeros(2,nuevasanot);
    anotaciondetectadaseg(2,:) = seg;
    if(anotextra)
        
        if(~isempty(anotaciondetectadaseg))
            
            anotaciondetectadaseg(2,end) = anotaciondetectadaseg(2,end) + 1;
       
        else
            
            anotaciondetectada(2,end) = anotaciondetectada(2,end) + 1;
            
        end
        
    end
    anotaciondetectada = [anotaciondetectada anotaciondetectadaseg];
    
    %Cantidad de nuevos latidos para añadirlos al segmento de comprobacion
    inilatseg = size(latidocorrecto,2);
    nuevoslat = size(QRSseg,1);
    latidocorrectoseg = zeros(2,nuevoslat);
    latidocorrectoseg(2,:) = seg;
    latidocorrecto = [latidocorrecto latidocorrectoseg];

    %Calcular el error de posicion de cada latido
    errorlatseg = -ones(2,nuevoslat);
    errorlatseg(2,:) = seg;
    errorlat = [errorlat errorlatseg];
    
    %Si hemos sustituido el ultima latido del segmento anterior
    %quitamos la posicion de la lista de latidos detectados
    if(sustultlat)
        posQRS(end) = [];
    end
    
    %Guardamos todos los latidos detectados
    posQRS = [posQRS QRSseg(:,1)']; 
    
    %Comprobamos si hay una anotacion de latido dentro de la ventana de
    %verificación de cada latido detectado
    for i=1:1:size(QRSseg,1)
        
        for j=1:1:length(lattime) 

            if(QRSseg(i,1) >= max(lattime(j)-mediavent,1) && (QRSseg(i,1) <= lattime(j) + mediavent))

                %Si se sutituye el ultimo latido del segmento anterior
                if(sustultlat && i==1)
                    %Marcamos el ultimo latido del segmento anterior como descartado
                    latidocorrecto(1,inilatseg) = 2;
                    %Quitamos la posicion de la lista de TP
                    posTP(end) = [];
                    %Volvemos a poner a cero la anotacion que se detectó 
                    %con el latido descartdo
                    if(ultlatdetsig)
                        anotaciondetectada(1,indultlatdet) = 0;
                    end
                    %Para que haga la comprobacion más rápido en las
                    %siguientes iteraciones
                    sustultlat = false;
                end
                
                latidocorrecto(1,inilatseg + i) = 1;
                posTP = [posTP QRSseg(i,1)];
                anotaciondetectada(1,inianotseg + j) = 1;
                
                %Calculamos el error en la detección
                errorlat(1,inilatseg + i) = abs(QRSseg(i,1) - lattime(j));
                
                %Guardamos el indice de la anotacion detectada con el
                %ultimo latido por si se descarta en el segmento siguiente
                if(i == size(QRSseg,1))
                    ultlatdet = true;
                    indultlatdet = inianotseg + j;
                end

                %Salimos del bucle de anotaciones
                break;

            end

        end
        
    end

    %Por si se descarta luego el ultimo latido y hay que quitar
    %tambien la posible anotacion detectada con él
    ultlatdetsig = ultlatdet;
    ultlatdet = false;
    
    %Comprobamos que ninguno de los candidatos que se expanden ya fue
    %seleccionado como el ultimo latido
    if(~isempty(QRScandsig) && ~isempty(QRSseg))
        elim = find((QRScandsig(:,1)+tamseg+tini) == QRSseg(end,1));
        QRScandsig(elim,:)=[];
    end
    
    if(~isempty(QRScandaltsig) && ~isempty(QRSseg))
        elim = find(QRScandaltsig(:,1)+tamseg+tini == QRSseg(end,1));
        QRScandaltsig(elim,:)=[];
    end
    
    %Candidatos que se expanden al segmento siguiente porque no están
    %enteros en éste
    QRScandsigant =  QRScandsig;
    QRScandaltsigant = QRScandaltsig;
    
end

%Volvemos al directorio original
cd(currentfolder);

%Resultados SIN tener en cuenta el estandar
%Si la ultima anotacion es del segmento posterior no procesado y no hay
%que tenerla en cuenta
if(ultimaanotacion == true)
    anotaciondetectada(:,size(anotaciondetectada,2))=[];
end
%Calcular los falsos negativos totales  
FNb = length(find(anotaciondetectada(1,:) == 0));
%Calcular los falsos positivos totales
FPb = length(find(latidocorrecto(1,:) == 0));
%Calcular los verdaderos positivos
TPb = length(find(anotaciondetectada(1,:) == 1));

%Calcular la sensibilidad total
Seb = TPb / (TPb + FNb) * 100;
%Calcular la predictividad positiva total
Ppb = TPb / (TPb + FPb) * 100;

%RMSE de las posiciones de los latidos
RMSE = sqrt((sum(errorlat.^2))/length(errorlat));

%Resultados siguiendo el estandar ANSI/AAMI
%Metadatos para el fichero de anotaciones
chanQRS = zeros(length(posQRS),1);
numQRS = zeros(length(posQRS),1);
%Todo 'N's
typeQRS = char(ones(length(posQRS),1) * 78);
%Todo '0's
subtypeQRS = char(ones(length(posQRS),1) * 48);
%Todo celdas vacías
commentQRS = cell(length(posQRS),1);

% NOTE: The WFDB Toolbox uses 0 based index, and MATLAB uses 1 based index.
%       Due to this difference annotation values ('ann') are shifted inside
%       this function in order to be compatible with the WFDB native
%       library. The MATLAB user should leave the indexing conversion to
%       the WFDB Toolbox.
senal = resultado(1).senal;
wrann([rutabd '\' senal],sufijo,posQRS',typeQRS,subtypeQRS,chanQRS,numQRS,commentQRS);

%Esta parte no se ejecuta siempre correctamente, en cocnreto la funcion bxb
%se queda colgada, por lo que se usa un procedimiento alternativo basado en
%las funciones originales de la wfdb ejecutadas mediante cygwin
%{
if(tinicomp == 0)
    tinicial = '0';
else
    tinicial = '300';
end
%Delete any previous temp.txt
if(exist('temp.txt', 'file'))
    delete('temp.txt');
end
tfinal = num2str(length(resultado(1).entrada)/frec);
resulcomp = bxb([rutabd '\' senal], 'atr', 'qrsd', 'temp.txt', tinicial, tfinal);

FN = sum(resulcomp.data(1:5,6));
FP = sum(resulcomp.data(6,1:5));
TP = sum(resulcomp.data(1:5,1));
Se = TP / (TP + FN) * 100;
Pp = TP / (TP + FP) * 100;
%}
