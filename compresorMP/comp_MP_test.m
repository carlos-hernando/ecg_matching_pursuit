function [CR,prd,prd1,prdseg,prd1seg,CRseg,RMSE,RMSEseg,maxidif,entrada,salida,Qfinal,coeficientes,atomos,codifseg,normascuant,mediascuant]=comp_MP_test...
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
%       distancia derecha; si convoluciones: centro, id semiatomo izquierdo,
%       id smiatomo derecho.
%
%    Parámetros de salida:
%       
%       CR: tasa de compresión de los diferenets metodos de compresion
%       prd: PRD de la señal reconstruida
%       prd1: otra medida de calidad eliminando la media
%       prdcc: medida de calidad con la linea base
%       prdseg: calidad PRD de los segmentos
%       prd1seg: calidad PRD1 de los segmentos
%       CRseg: compresion de los segmentos con los tres codificadores
%       RMSE: error cuadratico medio en mV
%       RMSEseg: error cuadratico medio de los segmentos
%       maxidif: error maximo en la recuperacion en mV
%       entrada: señal ECG original
%       salida: señal ECG recuperada
%       Qfinal: Valores finales de cuantificacion
%       coeficientes: coeficientes obtenidos con el MP (cada segmento es
%                     separado con un 0).
%       atomos: indices de los atomos [fichero, indice en el fichero] obtenidos 
%               con el MP (cada segmento es separado con un [0, 0]).
%       codifseg: si se usa la codificacion directa (0) o el codificador
%       entropico (1) en cada segmento
%       normascuant: normas de los segmentos cuantificadas
%       mediascuant: medias de los segmentos cuantificadas
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 31 mayo 2016           |
%             | AUTOR: Carlos Hernando Ramiro                |
%             |______________________________________________|
%
%

tic;

%Aqui se pone algo más bajo que el final de cuantificacion
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
ganancia =cabecera.gain(no_der);

%Guardamos el directorio actual para recuperarlo al final
currentfolder = pwd;

%Calculamos el tamaño de segmento según el diccionario empleado:

%Obtener los todos los ficheros de la carpeta diccio
ficheros=dir(diccio);
cd(diccio);

%Los dos primeros ficheros son . y ..
feval('load',ficheros(3).name);
%Al cargar el primer fichero del diccio
iant=3;
%Buscar el identificaor mas alto para despues codificar
maxparam = max(diccio.param(:));

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
coeficientes=[];
atomos=[];

%Reservar memoria
normas = zeros(1,num_segm);
normascuant = zeros(1,num_segm);
medias = zeros(1,num_segm);
mediascuant = zeros(1,num_segm);
Qfinal = zeros(1,num_segm);
prdseg = zeros(1,num_segm);
prd1seg = zeros(1,num_segm);
RMSEseg = zeros(1,num_segm);

for k=1:num_segm
  
  %Para depurar
  %numerosegmento = k;
    
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
          
         if i ~= iant
            feval('load',ficheros(i).name);
            
            %Buscar el identificador mas alto para despues codificar
            %solo con el primer segmento
            if(k==1)
                if(max(diccio.param(:))>maxparam)
                    maxparam = max(diccio.param(:));
                end
            end
            
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
         
         iant=i;
      
      end
         
      sdesnormed = sdesnormed - prodmax(iteracion)*prodatomo;
      
      recs =recs+(prodmax(iteracion)*prodatomo);
            
      recsfin=normascuant(k)*(recs+mediascuant(k));

      if tipo_cal==0
          calidad=prdecg(s,recsfin); %PRD
          calidad2=sqrt(sum(((recsfin-s)/ganancia).^2)/length(s)); %RMSE
          hacer = (iteracion<(segmento/2) && (abs(calidad-calidadobj)/calidadobj)>margencalidad && calidad>calidadobj);
      elseif tipo_cal==1
          calidad=sqrt(sum(((recsfin-s)/ganancia).^2)/length(s)); %RMSE
          calidad2=prdecg(s,recsfin); %PRD
          hacer = (iteracion<(segmento/2) && (abs(calidad-calidadobj)/calidadobj)>margencalidad && calidad>calidadobj);
      else
          calidad=prdecg(s,recsfin); %PRD
          calidad2=sqrt(sum(((recsfin-s)/ganancia).^2)/length(s)); %RMSE
          hacer = (iteracion < calidadobj);
      end

      %Para depurar:
      %calidadit(iteracion)=calidad;
      %calidad2it(iteracion)=calidad2;
      
  end
  
  clear diccio;
  iant = 0;
  
  %Cuantificamos los coeficientes
  [coefcuant, coefrec, Amax, Amin, recscuant, calidadseg, Qfin]=cuant_MP(ficheros, s, prodmax, Q, tipo_cal, calidadobj, mediascuant(k), normascuant(k), prodatomoind, precision);
  
  coeficientes=[coeficientes coefrec 0];
  atomos=[atomos; prodatomoind; 0 0 0];
  
  Qfinal(k)= Qfin;
  
  if (tipo_cal==0 || tipo_cal==2)
      prdseg(k)=calidadseg(1);
      prd1seg(k)=calidadseg(2);
      RMSEseg(k)=sqrt(sum(((recscuant-s)/ganancia).^2)/length(s)); %RMSE en mV
  elseif tipo_cal==1
      RMSEseg(k)=calidadseg(1);
      prdseg(k)=prdecg(s,recscuant); %PRD
      prd1seg(k)=calidadseg(2); %PRD1
  end
   
  %Bits de cabecera que hay que transmitir
  %16->norma, 16->media, 16->max cuantificacion, 
  %16->Min cuantificacion, 4->Q
  bitscabfij=16+16+16+16+4;
    
  %Metodo 1
  %Bits semi-fijos:
  %log2(segmento/2)->numero de atomos (como mucho la mitad del numero de
  %muestras menos uno)-> si todo 1's, entonces se va añadiendo otro grupo
  %de bits de tamaño log2(segmento/4) etc.
  numatom = length(coefcuant);
  indexbits = 1;
  newbitscab(1,k) = 0;
  
  while(numatom>0)
    newbitscab(1,k) = newbitscab(1,k) + ceil(log2(segmento/(2^indexbits)));
    numatom = numatom - ceil(segmento/(2^indexbits));
    indexbits = indexbits + 1;
  end
    
  bitscab(1,k) = bitscabfij + newbitscab(1,k);
  
  %Los coeficientes y los parametros de los atomos se codifican
  %directamente con el numero de bits necesarios
  bitscoef(1,k)=length(coefcuant)*Qfin;
  bitslado(1,k)=2*size(prodparam,1)*ceil(log2(maxparam));
  bitsdesp(1,k)=size(prodparam,1)*log2(segmento);
  
  codifseg(1,k) = 0;
  bits(1,k)=bitscab(1,k)+bitscoef(1,k)+bitslado(1,k)+bitsdesp(1,k);
  
  CRseg(1,k)=segmento*precision/bits(1,k);
  
  %Metodo 2
  
  %Los parametros de los atomos se codifican de manera aritmetica
  inputdatatest=cell(1,1);
  inputdatatest{1}=[prodparam(:,1) prodparam(:,2)];
  [codetest, compresultstest]=Arith06(inputdatatest);
  bitslado(2,k) = compresultstest(2,3);
  
  %Bits semi-fijos:
  %log2(segmento/2)->numero de bytes de la salida aritmetica (como mucho la mitad del numero de
  %muestras menos uno)-> si todo 1's, entonces se va añadiendo otro grupo
  %de bits de tamaño log2(segmento/4) etc.
  numatom = compresultstest(2,3)/8;
  indexbits = 1;
  newbitscab(2,k) = 0;
  
  while(numatom>0)
    newbitscab(2,k) = newbitscab(2,k) + ceil(log2(segmento/(2^indexbits)));
    numatom = numatom - ceil(segmento/(2^indexbits));
    indexbits = indexbits + 1;
  end
    
  bitscab(2,k) = bitscabfij + newbitscab(2,k);
  
  if((bitscab(2,k)+bitslado(2,k))<(bitscab(1,k)+bitslado(1,k)))
      codifseg(2,k) = 1;
      bits(2,k)=bitscab(2,k)+bitscoef(1,k)+bitslado(2,k)+bitsdesp(1,k);
  else
      codifseg(2,k) = 0;
      bits(2,k)=bitscab(1,k)+bitscoef(1,k)+bitslado(1,k)+bitsdesp(1,k);
  end
      
  %Bit extra que indica el tipo de codificacion
  bits(2,k) = bits(2,k) + 1;
      
  
  CRseg(2,k)=segmento*precision/bits(2,k);
  
  %Metodo 3
  
  %Los parametros de los atomos se codifican de manera Huffman
  [codetest, compresultstest]=Huff06(inputdatatest);
  bitslado(3,k) = compresultstest(2,3);
  
  %Bits semi-fijos:
  %log2(segmento/2)->numero de bytes de la salida aritmetica (como mucho la mitad del numero de
  %muestras menos uno)-> si todo 1's, entonces se va añadiendo otro grupo
  %de bits de tamaño log2(segmento/4) etc.
  numatom = compresultstest(2,3)/8;
  indexbits = 1;
  newbitscab(3,k) = 0;
  
  while(numatom>0)
    newbitscab(3,k) = newbitscab(3,k) + ceil(log2(segmento/(2^indexbits)));
    numatom = numatom - ceil(segmento/(2^indexbits));
    indexbits = indexbits + 1;
  end
    
  bitscab(3,k) = bitscabfij + newbitscab(3,k);
  
  if((bitscab(3,k)+bitslado(3,k))<(bitscab(1,k)+bitslado(1,k)))
      codifseg(3,k) = 1;
      bits(3,k)=bitscab(3,k)+bitscoef(1,k)+bitslado(3,k)+bitsdesp(1,k);
  else
      codifseg(3,k) = 0;
      bits(3,k)=bitscab(1,k)+bitscoef(1,k)+bitslado(1,k)+bitsdesp(1,k);
  end
  
  bits(3,k) = bits(3,k) + 1;
  
  CRseg(3,k)=segmento*precision/bits(3,k);
  
  %{
  %Pruebas en las que se pasan todos los parametros a bits, 
  %se reagrupan en elementos de 8 bits y se codifican
  bitstot = [];
  for n=1:length(coefcuant)
      bitstot = [bitstot de2bi(coefcuant(n),Qfin,'left-msb')];
  end
    
  for n=1:size(prodparam)
      bitstot = [bitstot de2bi(prodparam(n,1)-1,log2(max(max(prodparam(:,1:2)))),'left-msb')];
      bitstot = [bitstot de2bi(prodparam(n,2)-1,log2(max(max(prodparam(:,1:2)))),'left-msb')];
      bitstot = [bitstot de2bi(prodparam(n,3)-1,log2(segmento),'left-msb')];
  end
  
  aux = mod(length(bitstot),8);
 
  for n=1:aux
      bitstot = [bitstot 0];
  end
      
  for n=1:length(bitstot)/8
      findat(n) = bi2de(bitstot((n-1)*8+1:n*8));
  end
  
  for n=1:length(bitstot)/4
      findatb(n) = bi2de(bitstot((n-1)*4+1:n*4));
  end 
  
  datosent = cell(1,1);
  datosent{1} = findat(:);
  [code, compresults]=Arith06(datosent);
  recover = Arith06(code);
  bits(2,k)= bitscab+compresults(2,3);
  
  CRseg(2,k)=segmento*precision/bits(2,k);
  
  datosentb = cell(1,1);
  datosentb{1} = findatb(:);
  [codeb, compresultsb]=Arith06(datosentb);
  recoverb = Arith06(codeb);
  bits(3,k)= bitscab+compresultsb(2,3);
  
  CRseg(3,k)=segmento*precision/bits(3,k);
  
  [codec, compresultsc]=Huff06(datosent);
  recoverc = Huff06(codec);
  bits(4,k)= bitscab+compresultsc(2,3);
  
  CRseg(4,k)=segmento*precision/bits(4,k);
  
  [coded, compresultsd]=Huff06(datosentb);
  recoverd = Huff06(coded);
  bits(5,k)= bitscab+compresultsd(2,3);
  
  CRseg(5,k)=segmento*precision/bits(5,k);
  
  [codee, compresultse]=Arith07(datosent);
  recovere = Arith07(codee);
  bits(6,k)= bitscab+compresultse(2,3);
  
  CRseg(6,k)=segmento*precision/bits(6,k);
  
  [codef, compresultsf]=Arith07(datosentb);
  recoverf = Arith07(codef);
  bits(7,k)= bitscab+compresultsf(2,3);
  
  CRseg(7,k)=segmento*precision/bits(7,k);
  
  
  %Pruebas para ver si compensa agrupar los parametros de los atomos
  inputdatatesta=cell(1,1);
  %inputdatatesta{1}=prodparam(:);
  inputdatatesta{1}=[prodparam(:,1) prodparam(:,2)]
  [codetesta, compresultstesta]=Arith06(inputdatatesta);
  
  inputdatatestb=cell(1,2);
  %inputdatatestb{1}=[prodparam(:,1) prodparam(:,2)];
  %inputdatatestb{2}=prodparam(:,3);
  inputdatatestb{1}=prodparam(:,1);
  inputdatatestb{2}=prodparam(:,2);
  [codetestb, compresultstestb]=Arith06(inputdatatestb);
  recoverb = Arith06(codetestb);
  
  inputdatatestc=cell(1,1);
  inputdatatestd=cell(1,1);
  %inputdatatestc{1}=[prodparam(:,1) prodparam(:,2)];
  %inputdatatestd{1}=prodparam(:,3);
  inputdatatestc{1}=prodparam(:,1);
  inputdatatestd{1}=prodparam(:,2);
  [codetestc, compresultstestc]=Arith06(inputdatatestc);
  [codetestd, compresultstestd]=Arith06(inputdatatestd);
    
  juntos =  compresultstesta(2,3)
  semijuntos = compresultstestb(3,3)
  separados = compresultstestc(2,3) + compresultstestd(2,3) 
  %}
  
  signrec=[signrec recscuant];
  
end

CR(1) = length(signorg)*precision/sum(bits(1,:));
CR(2) = length(signorg)*precision/sum(bits(2,:));
CR(3) = length(signorg)*precision/sum(bits(3,:));
   
entrada=signorg;
salida=signrec;

prd=prdecg(entrada,salida);
prd1=prdmedio(entrada,salida);
entradamV=entrada/ganancia;
salidamV=salida/ganancia;
errormV=salidamV-entradamV;
RMSE=sqrt(sum(errormV.^2)/length(errormV));
maxidif=maxerror(entradamV, salidamV);

cd(currentfolder);
toc;
    
