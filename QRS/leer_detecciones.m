function [FN, FP, TP, Se, Pp, indi] = leer_detecciones(carpeta, sufijo)

%Guardamos el directorio actual para recuperarlo al final
diractual = pwd;

ficheros=dir(carpeta);

iteraciones = length(ficheros);

h = 0;

for i=3:iteraciones
    
    %Comprobamos que es un fichero de deteccion con el sufijo indciado
    if ~isempty(strfind(ficheros(i).name, '.txt')) && ~isempty(strfind(ficheros(i).name, 'det')) && ~isempty(strfind(ficheros(i).name, [sufijo '.']))
        
        tabla = importdata([carpeta '\' ficheros(i).name], ' ', 7);
        
        ref = strfind(ficheros(i).name, 'det');
        senal = str2num(ficheros(i).name(ref+3:ref+5));
        
        h = h + 1;
        
        indi(h,1) = senal;
        
        FN = sum(tabla.data(1:5,6));
        indi(h,2) = FN;
        FP = sum(tabla.data(6,1:5));
        indi(h,3) = FP;
        TP = sum(tabla.data(1:5,1));
        indi(h,4) = TP;
        Se = TP / (TP + FN) * 100;
        indi(h,5) = Se;
        Pp = TP / (TP + FP) * 100;
        indi(h,6) = Pp;
        
    end
    
end

FN = sum(indi(:,2));
FP = sum(indi(:,3));
TP = sum(indi(:,4));
Se = TP / (TP + FN) * 100;
Pp = TP / (TP + FP) * 100;

%Volvemos al directorio original
cd(diractual);