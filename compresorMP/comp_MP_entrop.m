function [CRlow,CRhigh,entroplow,entrophigh,prd,prd1,prdcc,prdseg,prd1seg,CRlowseg,CRhighseg,RMSE,RMSEseg,maxidif,entrada,salida,entradamV,salidamV,Qfinal,coeficientes,atomos]=comp_MP_entrop...
(senal,no_der,long_signal,tipo_cal,calidadobj,Q,ultimo,diccio)

% Función que realiza el calculo de compresión maxima, con la entropía, de una señal ECG  
% usando un diccionario de atomos y la tecnica Matching Pursuit con
% convoluciones circulares mediante FFT
%
%    Parámetros de entrada:
%       senal: string de caracteres indicando el nombre de la señal a procesar.
%       no_der: entero que indica cual de las dos derivaciones del fichero del MIT
%               se va a emplear.
%       long_signal: longitud que se va a procesar de la señal en segundos.
%       tipo_cal: '0'->PRD, '1'-> RMSE, '2'->iteraciones
%       calidadobj: calidad con que se quiere recuperar la señal. PRD en %
%       o RMSE en mV o numero de iteraciones.
%       Q: bits de cuantificación de los coeficientes del MP, si Q==0
%       Q adaptativo, menos con iteracioens que Q=9
%       ultimo: si se tiene en cuenta (1) o no (0) el ultimo segmento
%       diccio: carpeta con los ficheros del diccionario a aplicar, son del tipo
%       atom_dic: atomos a utilizar, su longitud debe ser
%       igual al segmento y param_dic: parametros que representan los atomos 
%       del diccionario: si triangulos: centro, distancia izquierda, 
%       distancia derecha; si convoluciones: centro, long inicial izq, 
%       convoluciones izq, long ini derecha y conv derecha.
%
%    Parámetros de salida:
%       
%       CRlow: tasa de compresión simulada de manera realista
%       CRhigh: tasa de compresion simulada de manera optimista
%       prd: PRD de la señal reconstruida
%       prd1: otra medida de calidad eliminando la media
%       prdcc: medida de calidad con la linea base
%       prdseg: calidad PRD de los segmentos
%       prd1seg: calidad PRD1 de los segmentos
%       CRlowseg: compresion de los segmentos realista
%       CRhighseg: compresion de los segmentos optimista
%       RMSE: error cuadratico medio en mV
%       RMSEseg: error cuadratico medio de los segmentos
%       maxidif: error maximo en la recuperacion en mV
%       entrada: señal ECG original
%       salida: señal ECG recuperada
%       entradamV: señal ECG original en mV
%       salidamV: señal ECG recuperada en mV
%       Qfinal: Valores finales de cuantificacion
%       coeficientes: coeficientes obtenidos con el MP (cada segmento es
%                     separado con un 0).
%       atomos: indices de los atomos [fichero, indice en el fichero] obtenidos 
%               con el MP (cada segmento es separado con un [0, 0]).
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 9 febrero 2015         |
%             | AUTOR: Carlos Hernando Ramiro                |
%             |______________________________________________|
%
%

tic;

margencalidad=0.05;

y1=[senal '.hea'];
y2=[senal '.dat'];
cabecera=feval('readheader_nsrdb',y1);
fid=feval('fopen',y2,'r');
[datos,cerrado]=lee_MIT_segundos(fid,long_signal,cabecera,no_der);
signorg=datos-cabecera.adczero(no_der);% A la señal original se le ajusta el cero.
if cerrado
    fclose(fid);
end

precision=cabecera.adcres(no_der);
fs=cabecera.freq;

%Guardamos el directorio actual para recuperarlo al final
currentfolder = pwd;

%Calculamos el tamaño de segmento según el diccionario empleado:

%Obtener los todos los ficheros de la carpeta diccio
ficheros=dir(diccio);
cd(diccio);

%Los dos primeros ficheros son . y ..
feval('load',ficheros(3).name);

segmento=size(diccio.atomos,1);

num_segm=ceil((length(signorg)/segmento));
n_ceros_mas=(num_segm*segmento)-length(signorg);

signorg=[signorg zeros(1,n_ceros_mas)];

%En este caso NO tenemos en cuenta el ultimo segmento
if n_ceros_mas>0 && ultimo==0
    signorg=signorg(1:length(signorg)-segmento);
    num_segm=num_segm-1;
end

signrec=[];
Qfinal=[];

coeficientes=[];
atomos=[];

for k=1:num_segm
  
  numerosegmento = k
    
  prodparam=[];
  prodmax=[];
  prodatomo=[];
  prodatomoind=[];
  iteracion=0;
  recs=zeros(1,segmento);
  
  s=signorg((segmento*(k-1)+1):segmento*(k));
  
  normas(k) = norm(s);
  
  %Cuantificamos la norma con 16 bits
  Mmquantizer = quantizer([16,5],'float');
  normascuant(k) = quantize(Mmquantizer,normas(k));    
  
  sdesnor=s/normascuant(k);
  
  medias(k)=mean(sdesnor);
  
  %Cuantificacmos la media con 16 bits
  mediascuant(k) = quantize(Mmquantizer,medias(k));    
  
  sdesnormed=sdesnor-mediascuant(k);
  
  hacer = true;
  while hacer

      iteracion=iteracion+1;
      prodmax(iteracion)=0;
      
      %Los dos primeros ficheros son . y ..
      for i=3:length(ficheros) 
          
         if iteracion == 1 || (i ~= 3)
            feval('load',ficheros(i).name);
         end
         
         %Meto directamente aux1 y aux2 en bsxfun y luego
         %aux3 en ifft para que no se guarden en memoria
         %Hace los mismo que cconv pero para una matriz
         %aux1=fft(sdesnormed,segmento).';
         %aux2=fft(diccio.flipatomos,segmento);
         %Multiplicar el vector por todas las columnas
         %aux3=bsxfun(@times,fft(sdesnormed,segmento).',fft(diccio.flipatomos,segmento));
         
         prods = ifft(bsxfun(@times,fft(sdesnormed,segmento).',fft(diccio.flipatomos,segmento)));
         [prodmaxaux, index] = max(abs(prods(:)));
                      
         if prodmaxaux>abs(prodmax(iteracion))

             [jprod,iprod]=ind2sub(size(prods),index(1));
             prodmax(iteracion) = prods(index);
             %Añadimos como ultimo parametro la muestra donde comienza 
             %el atomo [1-256]
             prodparam(iteracion,:)=[diccio.param(:,iprod)' jprod];
             %prodsenal(iteracion,:)=diccio.atomos(jprod(1),:); %para
             %ver el atomo y depurar
             prodatomo=circshift(diccio.atomos(:,iprod),(jprod-1))';
             %Almacenamos: numero de fichero, numero de atomo y
             %desplazamiento
             prodatomoind(iteracion,:)=[i iprod (jprod-1)];

         end
         
         clear prods;
         %clear diccio;
      
      end
         
      sdesnormed = sdesnormed - prodmax(iteracion)*prodatomo;
      
      recs =recs+(prodmax(iteracion)*prodatomo);
            
      recsfin=normascuant(k)*(recs+mediascuant(k));
      
      if tipo_cal==0
          calidad=prdecg(s,recsfin); %PRD
          calidad2=sqrt(sum(((recsfin-s)*10/(2^precision-1)).^2)/length(s)); %RMSE
          hacer = (iteracion<(segmento/2) && (abs(calidad-calidadobj)/calidadobj)>margencalidad && calidad>calidadobj);
      elseif tipo_cal==1
          calidad=sqrt(sum(((recsfin-s)*10/(2^precision-1)).^2)/length(s)); %RMSE
          calidad2=prdecg(s,recsfin); %PRD
          hacer = (iteracion<(segmento/2) && (abs(calidad-calidadobj)/calidadobj)>margencalidad && calidad>calidadobj);
      else
          calidad=prdecg(s,recsfin); %PRD
          calidad2=sqrt(sum(((recsfin-s)*10/(2^precision-1)).^2)/length(s)); %RMSE
          hacer = (iteracion < calidadobj);
      end
      
      calidadit(iteracion)=calidad;
      calidad2it(iteracion)=calidad2;
      
  end
  
  clear diccio;
  
  %Cuantificamos los coeficientes
  [coefcuant, coefrec, Amax, Amin, recscuant, calidadseg, Qfin]=cuant_MP(ficheros, s, prodmax, Q, tipo_cal, calidadobj, mediascuant(k), normascuant(k), prodatomoind, precision);
  
  coeficientes=[coeficientes coefrec 0];
  atomos=[atomos; prodatomoind; 0 0 0];
  
  Qfinal=[Qfinal Qfin];
  
  if (tipo_cal==0 || tipo_cal==2)
      prdseg(k)=calidadseg(1);
      prd1seg(k)=calidadseg(2);
      RMSEseg(k)=sqrt(sum(((recscuant-s)*10/(2^precision-1)).^2)/length(s)); %RMSE en mV
  elseif tipo_cal==1
      RMSEseg(k)=calidadseg;
      prdseg(k)=prdecg(s,recscuant); %PRD
      prd1seg(k)=prdmedio(s,recscuant); %PRD1
  end
   
  %CR calculado de manera realista
  %16->norma, 16->media, 16->max cuantificacion, 16->Min cuantificacion
  %4->Q, (length(coefcuant)*Qfin)->coeficientes cuantificados
  bitslow(k)=16+16+16+16+4+(length(coefcuant)*Qfin);
 
  %Bits de los parametros de los atomos      
  [prob,entropparam(k),simbolos]=calc_entrop(prodparam(:));
  bitsparamlow(k)=ceil(entropparam(k)*length(prodparam(:)));
  bitslow(k)=bitslow(k)+bitsparamlow(k);
  
  entroplow(k,1)=Qfin;
  entroplow(k,2)=entropparam(k);
  
  CRlowseg(k)=segmento*precision/bitslow(k);
  
  
  %CR calculado de manera optimista
  %Bits necesarios para descomprimir
  bitshigh(k)=16+16+16+16+4; %Norma, media, Max cuantificacion, Min cuantificacion, Q
  
  %Bits de los coeficientes cuantificados
  [prob2,entropcoef2(k),simbolos2]=calc_entrop(coefcuant);
  bitscoefhigh(k)=ceil(entropcoef2(k)*length(coefcuant));
  bitshigh(k)=bitshigh(k)+bitscoefhigh(k);
 
  %Bits de los parametros de los atomos
  for i=1:size(prodparam,2)
  
    [prob2,entropparam2(k,i),simbolos2]=calc_entrop(prodparam(:,i));
    bitsparamhigh(k,i)=ceil(entropparam2(k,i)*length(prodparam(:,i)));
    bitshigh(k)=bitshigh(k)+bitsparamhigh(k,i);
  
  end
  
  entrophigh(k,1)=entropcoef2(k);
  entrophigh(k,2)=entropparam2(k,1);
  entrophigh(k,3)=entropparam2(k,2);
  entrophigh(k,4)=entropparam2(k,3);
  
  CRhighseg(k)=segmento*precision/bitshigh(k);
  
  signrec=[signrec recscuant];
  
end

CRlow=length(signorg)*precision/sum(bitslow);
CRhigh=length(signorg)*precision/sum(bitshigh);

entrada=signorg;
salida=signrec;
    
prd=prdecg(entrada,salida);
prd1=prdmedio(entrada,salida);
prdcc=prdecg(entrada+cabecera.adczero(no_der),salida+cabecera.adczero(no_der));
entradamV=entrada*10/(2^precision-1);
salidamV=salida*10/(2^precision-1);
errormV=salidamV-entradamV;
RMSE=sqrt(sum(errormV.^2)/length(errormV));
maxidif=maxerror(entradamV, salidamV);

cd(currentfolder);
toc;
    
