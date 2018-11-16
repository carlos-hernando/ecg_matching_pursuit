function [datos,cerrado]=lee_MIT_segundos(fid,long_signal,cabecera,no_der);

% Funcion que lee un intervalo (en segundos) del formato 212 de los
% ficheros del MIT-BIH desde la posición en que se encuentre el puntero
% en ese momento.
%
%          [datos,cerrado]=lee_MIT_Segundos(fid,long_signal,cabecera,no_der);
%
%    Parámetros de entrada:
%              fid: indentificador del fichero.
%              long_signal: nº entero que designa la duración en segundos del tramo de señal.
%              cabecera: variable de cabecera que indica las características del fichero.
%              no_der: entero que designa la derivación que se quiere capturar.
%    Parámetros de salida:
%              datos: señal que ocupa el intervalo seleccionado.
%              cerrado: mientras esta variable sea igual a 2, el fichero permanecerá abierto.


if ((cabecera.nsamp)-(ftell(fid)/3))>=long_signal*cabecera.freq
   a=fread(fid,[3 long_signal*cabecera.freq],'uchar');
   cerrado=2;
else 
   a=fread(fid,[3 cabecera.nsamp-(ftell(fid)/3)],'uchar');
   cerrado=fclose(fid);
end


if no_der==1
   b1=bitand(a(2,:),15)*2^8;
   datos=a(1,:)+b1;
else
   b2=bitand(a(2,:),240)*2^4;
   datos=a(3,:)+b2;
end
   
neg=bitget(datos,12);% Se buscan aquellas muestras negativas.
datos(find(neg))=-1*(bitxor(datos(find(neg)),2^12-1)+1);% C2 de las muestras negativas