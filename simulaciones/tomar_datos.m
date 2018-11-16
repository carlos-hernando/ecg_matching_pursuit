function [prd, CR, CR1, CR2, CR3, RMSE, NSA, MAX, codifdir1, codifdir2, codifdir3]=tomar_datos(carpeta)

%Guardamos el directorio actual para recuperarlo al final
currentfolder = pwd;

ficheros=dir(carpeta);
cd(carpeta);

i=3;
h=1;

flag=0;

while i<=length(ficheros) && flag==0
    
    %Comprobamos que es un fichero .mat no .fig u otros
    if strfind(ficheros(i).name, '.mat')
        
        feval('load',ficheros(i).name); 
     
        for j=1:length(resultado)

            prd(h,j)= resultado(j).prd;
            RMSE(h,j)= resultado(j).RMSE;
            if size(resultado(j).CR,2)==3
                CR1(h,j)= resultado(j).CR(1);
                CR2(h,j)= resultado(j).CR(2);
                CR3(h,j)= resultado(j).CR(3);
                CR(h,j)=0;
                NSA(h,j)= length(resultado(j).entrada)/length(resultado(j).coeficientes);
                MAX(h,j)=resultado(j).maxidif;
                codifdir1(h,j)=length(resultado(j).codifseg(1,:))-sum(resultado(j).codifseg(1,:));
                codifdir2(h,j)=length(resultado(j).codifseg(2,:))-sum(resultado(j).codifseg(2,:));
                codifdir3(h,j)=length(resultado(j).codifseg(3,:))-sum(resultado(j).codifseg(3,:));
            else
                CR(h,j) = resultado(j).CR;
                CR1(h,j)= 0;
                CR2(h,j)= 0;
                CR3(h,j)= 0;
                NSA(h,j)=0;
                MAX(h,j)=resultado(j).maxidif*10/(2^11-1);
                codifdir1(h,j)=0;
                codifdir2(h,j)=0;
                codifdir3(h,j)=0;
            end
         
        end

        h=h+1;
    end
    
    i=i+1;

end

%Volvemos al directorio original
cd(currentfolder);