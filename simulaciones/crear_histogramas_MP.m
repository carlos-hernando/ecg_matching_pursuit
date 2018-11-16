function [prd, histdata]=crear_histogramas_MP(carpeta, diccio)

diractual = pwd;

ficheros=dir(carpeta);
cd(carpeta);

i=3;

while i<=length(ficheros)

    %Comprobamos que es un fichero .mat no .fig u otros
    if strfind(ficheros(i).name, '.mat')

        feval('load',ficheros(i).name); 

        for j=1:length(resultado)

            prdtod(i-2,j)= resultado(j).prd;
            atomostod{i-2,j} = resultado(j).atomos;
            coeficientestod{i-2,j}= resultado(j).coeficientes;

        end

    end

        i=i+1;
    
end

if size(prdtod,1)>1
    prd = mean(prdtod);
else
    prd = prdtod;
end

atomos=[];
coeficientes=[];

for i=1:size(prdtod,2)
    
    atomos=[];
    coeficientes=[];
    
    for j=1:size(prdtod,1)
        atomos = [atomos; atomostod{j,i}];
        coeficientes = [coeficientes coeficientestod{j,i}];
    end
    
    [histdata(i).longlados, histdata(i).val, histdata(i).valcoef]=obtener_longlados(atomos, coeficientes, diccio);
        
end

cd(diractual);
