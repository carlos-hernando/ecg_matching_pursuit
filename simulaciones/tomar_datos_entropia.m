function [prd, CRlow, CRhigh, RMSE, NSA, MAX]=tomar_datos_entropia(carpeta)

ficheros=dir(carpeta);
cd(carpeta);

i=3;
h=1;

while i<=length(ficheros)
    
    %Comprobamos que es un fichero .mat no .fig u otros
    if strfind(ficheros(i).name, '.mat')
        
        feval('load',ficheros(i).name); 
     
        for j=1:length(resultado)

            prd(h,j)= resultado(j).prd;
            CRlow(h,j)= resultado(j).CRlow;
            CRhigh(h,j)= resultado(j).CRhigh;
            RMSE(h,j)= resultado(j).RMSE;
            NSA(h,j)= length(resultado(j).entrada)/length(resultado(j).coeficientes);
            MAX(h,j)=resultado(j).maxidif;
         
        end

        h=h+1;
    end
    
    i=i+1;

end