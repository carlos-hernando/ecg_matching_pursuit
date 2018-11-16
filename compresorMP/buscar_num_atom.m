function []=buscar_num_atom(carpeta, segmento)
%Se obtienen los atomos por segmento de todos los ficheros de resultados de
%compresión MP de la carpeta de entrada


%Guardamos el directorio actual para recuperarlo al final
currentfolder = pwd;

ficheros=dir(carpeta);
cd(carpeta);

carpeta

i=3;
h=1;

flag=0;

while i<=length(ficheros) && flag==0
    
    %Comprobamos que es un fichero .mat no .fig u otros
    if strfind(ficheros(i).name, '.mat')
        
        feval('load',ficheros(i).name);
     
        for j=1:length(resultado)
            
            a = find(resultado(j).coeficientes==0);
            b = a - [0 a(1:end-1)];
            b = b - 1; 
            c = find(b >= segmento/2);           
            
            
            if(length(c)>0)
                errormsg = 'Error!!!'
                senal = resultado(j).senal
                calidad = resultado(j).calidadobj
                seg_error = c
                valores = b(c)                
            end
            
        end

        h=h+1;
        
    end
    
    i=i+1;

end

%Volvemos al directorio original
cd(currentfolder);