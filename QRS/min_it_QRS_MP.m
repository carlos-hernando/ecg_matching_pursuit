function [minit, indi] = min_it_QRS_MP(senales)
%Calcula el numero minimo de iteraciones obtenidas teniendo en cuenta todos
%los segmentos para cada una de la señales de la carpeta de entrada

%Guardamos el directorio actual para recuperarlo al final
diractual = pwd;

ficheros=dir(senales);
cd(senales);

h = 0;
i = 3;

while i<=length(ficheros)
    
    %Comprobamos que es un fichero .mat no .fig u otros
    if strfind(ficheros(i).name, '.mat')
        
        ficheros(i).name
        
        ficherosenal = [senales '\' ficheros(i).name]; 
     
        load(ficherosenal)
        
        inter = find(resultado(1,1).coeficientes == 0);
        interref = [0 inter(1:end-1)];
        itpersegment = inter - interref;
        itpersegment = itpersegment - 1;
               
        h = h + 1;
       
        indi(h) = min(itpersegment);
        
    end
    
    i = i + 1;
    
end

minit = min(indi);

%Volvemos al directorio original
cd(diractual);
        