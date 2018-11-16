function [coefcuant, coefrec, Amax, Amin, recscuant, calidadseg, Q]=cuant_MP(ficheros, s, coefentrada, Q, tipo_cal, calidadobj, mediascuant, normascuant, prodatomoind, ganancia)

% Función que realiza la cuantificacion de los coeficientes en el caso de 
% haberlos obtenido con el metodo de Matching Pursuit de convoluciones
% circulares
%
%    Parámetros de entrada:
%       ficheros: archivos de la carpeta del diccionario
%       s: segmento original
%       coefentrada: coeficientes obtenidos con el MP
%       Q: bits de cuantificación de los coeficientes del MP, si 
%       calidad por iteraciones y Q==0 se usan 9 bits. 
%       tipo_cal: '0'->PRD, '1'-> RMSE, '2'-> iteraciones
%       calidadobj: calidad con que se quiere recuperar la señal. PRD en %
%       o RMSE en mV o numero de iteraciones.
%       mediascuant: valor medio del segmento
%       normascuant: norma del segmento
%       prodatomoind: indices de los atomos usados en cada iteracion con el
%       formato: numero de fichero, fila del atomo y desplazamiento
%       ganancia: ganancia de adquisicion de la señal para calcular el RMSE
%
%    Parámetros de salida:
%       
%       coefcuant: coeficientes cuantificados
%       coefrec: coefcientes recuperados tars cuantificacion
%       Amax: valor maximo de cuantificacion
%       Amin: valor minimo de cuantificacion
%       recscuant: segmento recuperado
%       calidadseg: calidad obtenida en el segmento
%       Q: numero de bits de cuantificacion empleado
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 31 de mayo 2016        |
%             | AUTOR: Carlos Hernando Ramiro                |
%             |______________________________________________|
%
%

Mmquantizer = quantizer([16,5],'float');
Mmquantizer2=quantizer([16,5],'float','ceil');
Aminsin=min(coefentrada);
Amin=quantize(Mmquantizer, Aminsin);
Amaxsin=max(coefentrada);
Amax=quantize(Mmquantizer2, Amaxsin);

margencalidad=0.1;

%Si no se ha puesto Q fijo a pesar de usar iteraciones, se fija Q
if tipo_cal==2 && Q==0
    Q=9;
end

if Q==0
  
    %Minimo Q=6
    Q=5;
    
    hacer = true;
    while hacer
    
        Q=Q+1;
        
        %Límite maximo Q=11
        if (Q>11)
                break
        end
       
        coefcuant=round(((2^Q-1)*((coefentrada-Amin)/(Amax-Amin))));
        
        coefrec=((coefcuant*(Amax-Amin))/(2^Q-1))+Amin;
        
        recscuant=zeros(1,length(s));
        
        for i=1:length(coefentrada)
            if i==1 || (prodatomoind(i,1)~=prodatomoind(i-1,1))
                feval('load',ficheros(prodatomoind(i,1)).name);
            end
            prodatomo=circshift(diccio.atomos(:,prodatomoind(i,2)),prodatomoind(i,3))';
            recscuant = recscuant+(coefrec(i)*prodatomo);
        end
      
        recscuant=normascuant*(recscuant+mediascuant);
        
        if (tipo_cal==0 || tipo_cal==2)
            calidadseg=prdecg(s,recscuant); %PRD
        else
            calidadseg=sqrt(sum(((recscuant-s)/ganancia).^2)/length(s)); %RMSE
        end
        
     %Que sea mayor que calidadobj pero dentro del margen de error o que
     %sea menor que calidadobj
     hacer = ((abs(calidadseg-calidadobj)/calidadobj)>margencalidad && calidadseg>calidadobj);
        
    end
else
    
    coefcuant=round(((2^Q-1)*((coefentrada-Amin)/(Amax-Amin))));

    coefrec=((coefcuant*(Amax-Amin))/(2^Q-1))+Amin;
    
    recscuant=zeros(1,length(s));
        
    for i=1:length(coefentrada)
        if i==1 || (prodatomoind(i,1)~=prodatomoind(i-1,1))
            feval('load',ficheros(prodatomoind(i,1)).name);
        end
        prodatomo=circshift(diccio.atomos(prodatomoind(i,2),:)',prodatomoind(i,3))';
        recscuant = recscuant+(coefrec(i)*prodatomo);
    end
      
    recscuant=normascuant*(recscuant+mediascuant);
        
    if (tipo_cal==0 || tipo_cal==2)
        calidadseg=prdecg(s,recscuant); %PRD
    else
        calidadseg=sqrt(sum(((recscuant-s)/ganancia).^2)/length(s)); %RMSE
    end
          
end

calidadseg(2)=prdmedio(s,recscuant); %PRD1



