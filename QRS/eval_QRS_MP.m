function [indi, total] = eval_QRS_MP(senales, diccio, idcal, segini, segfin, itfin, frecuencia, ganancia, senalref, rutabd, tinicomp, sufijo)

%Guardamos el directorio actual para recuperarlo al final
diractual = pwd;

ficheros=dir(senales);

h = 0;
i = 3;

FNtotal = 0;
FPtotal = 0;
TPtotal = 0;

while i<=length(ficheros)
    
    %Comprobamos que es un fichero .mat no .fig u otros
    if strfind(ficheros(i).name, '.mat')
        
        ficheros(i).name
        
        ficherosenal = [senales '\' ficheros(i).name]; 
     
        [FN, FP, TP, Se, Pp] = detectar_QRS_MP(ficherosenal, diccio, idcal, segini, segfin, itfin, frecuencia, ganancia, senalref, rutabd, tinicomp, sufijo);
        
        %Obtener el id de la senal
        ref1 = strfind(ficheros(i).name,'senal');
        ref2 = strfind(ficheros(i).name,'tiempo');
        senaltexto = ficheros(i).name(ref1+6:ref2-2); 
        
        h = h + 1;
        
        indi(h, 1) = str2num(senaltexto);
        indi(h, 2) = FN;
        indi(h, 3) = FP;
        indi(h, 4) = TP;
        indi(h, 5) = Se;
        indi(h, 6) = Pp;
        
        FNtotal = FNtotal + FN;
        FPtotal = FPtotal + FP;
        TPtotal = TPtotal + TP;
        
    end
    
    i = i + 1;
    
end

total = struct('diccio',[],'latidos',[],'TP',[],'FN',[],'FP',[],...
                'Se',[],'Ppos',[],'media',[],'minimo',[]);

%Obtener el nombre del diccionario y el resto de detos
ref = strfind(diccio,'dic_');
total.diccio = diccio(ref+4:length(diccio));
total.latidos = TPtotal + FNtotal;
total.TP = TPtotal;
total.FN = FNtotal;
total.FP = FPtotal;
total.Se = TPtotal / (TPtotal + FNtotal) * 100;
total.Ppos = TPtotal / (TPtotal + FPtotal) * 100;
total.media = mean([total.Se total.Ppos]);
total.minimo = min([total.Se total.Ppos]);

%Volvemos al directorio original
cd(diractual);
        