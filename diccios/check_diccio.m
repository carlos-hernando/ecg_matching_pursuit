function[tamano, error]=check_diccio(dirdiccio)
%Funcion que comprueba el numero de atomos del diccionario y si
%hay atomos repetidos en los ficheros que lo componen (lo ideal seria
%comprobar si hay repetidos en general pero de momento no se como hacerlo)
%dirdiccio: direccion completa de la carpeta que contiene el diccionario

%Obtener los todos los ficheros de la carpeta diccio
ficheros=dir(dirdiccio);
cd(dirdiccio);

tamano = 0;
error = [];
division = 2^16;

for i=3:length(ficheros)
    
    i-2
    
    feval('load',ficheros(i).name);
    tamanoparte = size(diccio.atomos,1);
    tamano = tamano + tamanoparte;
    
    k=1;
    while (tamanoparte >0)
           
        if(tamanoparte >= division)
            partediccio = diccio.atomos(division*(k-1)+1:division*k,:);
            tamanoparte = tamanoparte - division;
        else
            partediccio = diccio.atomos(division*(k-1)+1:division*(k-1)+tamanoparte,:);
            tamanoparte = 0;
        end
            
            [uniquematrix,uniqueindeces] = unique(partediccio,'rows','first');  %Finds indices of unique rows
            repeatedindeces = setdiff(1:size(partediccio,1),uniqueindeces);  %Finds indices of repeats
            

        if ~isempty(repeatedindeces)
            if k==1
                error = [error i-2 repeatedindices];
            else
                error = [error repeatedindices];
            end
        end
        
        
    end    
            
   clear diccio;  
end