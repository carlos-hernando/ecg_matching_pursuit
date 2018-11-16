function resultados = eval_QRS_MP_todos(carpresul, carpdic, idcal, segini, segfin, itfin, frecuencia, ganancia, senalref, rutabd, tinicomp, sufijos)

%sufijos tiene que tener el formato:
% [{'qrs'} {'qrsd'}]

%Guardamos el directorio actual para recuperarlo al final
diractual = pwd;

%Estructura con los resultados
resultados = struct('idcal', {}, 'segini', {}, 'segfin', {},...
                'itfin', {}, 'frecuencia', {}, 'senalref', {}, 'resdic', {}, 'indi', {});

ficheros=dir(carpresul);

i = 3;
h = 0;
iteraciones = length(ficheros);

while i<=iteraciones
    
    %Obtener el nombre del diccionario
    ref = strfind(ficheros(i).name,'dic_');
    ref2 = strfind(ficheros(i).name,'tip');
    nomdic = ficheros(i).name(ref+4:ref2-2);
    diccio = [carpdic '\dic_' nomdic];
    diccios{i-2}=nomdic;
    
    iteracion = i - 2
    disp('de')
    iteraciones - 2
    nomdic
    disp('con sufijo')
    sufijos{i-2}
    
    [indi, total] = eval_QRS_MP([carpresul '\' ficheros(i).name], diccio, idcal, segini, segfin, itfin, frecuencia, ganancia, senalref, rutabd, tinicomp, sufijos{i-2});
     
    h = h + 1;
    resultados(h).idcal = idcal;
    resultados(h).segini = segini;
    resultados(h).segfin = segfin;
    resultados(h).itfin = itfin;
    resultados(h).frecuencia = frecuencia;
    resultados(h).senalref = senalref;
    resultados(h).resdic = total;
    resultados(h).indi = indi;
    
    i = i + 1;
     
end

for i=1:length(diccios)
    disp(diccios{i})
    disp(sufijos{i})
end

%Volvemos al directorio original
cd(diractual);
     
