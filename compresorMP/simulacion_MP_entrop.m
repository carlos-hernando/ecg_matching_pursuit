function[]=simulacion_MP_entrop(gruposenal,diccio,tipo_cal,calidad,Q)

% Funcion que simula el procesdo de señales con el compresor Matching Pursuit
% y guarda lso resulatdos para cada señal en un fichereo independiente
% 
%    Parámetros de entrada:
%       gruposenal: (0) señales completas de un fichero de texto determinado
%       (1) 11 señales de 2 minutos (2) 11 señales de 10 min (3) 19 señales
%       de 1 minuto (4) cualqueir otro caso particular
%       diccio: ruta completa de la carpeta que contiene el diccionario
%       tipo_cal: '0'->PRD, '1'-> RMSE, '2'->iteraciones
%       calidadobj: calidad con que se quiere recuperar la señal. PRD en %
%       o RMSE en mV o numero de iteraciones.
%       Q: bits de cuantificación de los coeficientes del MP, si Q==0
%       Q adaptativo, menos con iteracioens que Q=9.
%
%    Parámetros de salida: N/A
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 15 de enero de 2016    |
%             | AUTOR: Carlos Hernando                       |
%             |______________________________________________|
%


%Laboratorio
%cd('C:\Documents and Settings\Carlos B\Escritorio');
%Portatil
%cd('C:\Users\carlos\Desktop');

if (gruposenal ~= 1 & gruposenal ~= 2 & gruposenal ~= 3 & gruposenal ~= 4)
    grupo=0;
else
    grupo=gruposenal;
end

switch grupo
    case 0
        k=1;
        fid=feval('fopen',gruposenal,'r');
        tline{k} = fgetl(fid);
        while ischar(tline{k})
            k=k+1;
            tline{k} = fgetl(fid);
        end
        fclose(fid);
        
        for i=1:length(tline)-1
            ficheros{i}=tline{i};
            senales(i)=str2num(ficheros{i});
        end
        
        etiq=gruposenal;
        
    case 1
        senales=[100 101 102 103 107 109 111 115 117 118 119];
        tiempo=[120];
        etiq='GRUPOde100a119_2min';
        
    case 2
        senales=[100 101 102 103 107 109 111 115 117 118 119];
        tiempo=[600];
        etiq='GRUPOde100a119_10min';
        
    case 3
        senales=[104 107 111 112 115 116 117 118 119 201 207 208 209 212 213 214 228 231 232];
        tiempo=[60];
        etiq='GRUPOde104a232_1min';
        
    case 4
        senales=[117];
        tiempo=[120];
        etiq='SENAL117_2min';
        
end

%Obtenemos el nombre del diccionario
k = strfind(diccio, '\');
labeldiccio=diccio(k(length(k))+1:length(diccio));

resultado = struct('senal',{},'no_der',{},'long_signal',{},'tipo_cal',{},'calidadobj',{},...
    'Q',{},'ultimo',{},'diccio',{}, ...
    'CRlow',{},'CRhigh',{},'prd',{},'prd1',{},'prdcc',{},'prdseg',{},'prd1seg',{},'CRlowseg',{},'CRhighseg',{}, ...
    'RMSE',{},'RMSEseg',{},'maxidif',{},'entrada',{},'salida',{},'Qfinal',{},'coeficientes',{},'atomos',{});
  
trans='MP';

feval('mkdir',[trans '_' etiq '_' labeldiccio '_tip_' num2str(tipo_cal) '_cal_' num2str(floor(min(min(calidad)))) ...
    '_' num2str(ceil(max(max(calidad)))) '_Q_' num2str(Q)]);
oldfolder=feval('cd',[trans '_' etiq '_' labeldiccio '_tip_' num2str(tipo_cal) '_cal_' num2str(floor(min(min(calidad)))) ...
    '_' num2str(ceil(max(max(calidad)))) '_Q_' num2str(Q)]);

for i=1:length(senales)
    
    ticgrupo=tic;
    
    if grupo==0
        cabecera=feval('readheader_nsrdb',[ficheros{i} '.hea']);
        tiempo=cabecera.nsamp/cabecera.freq;
    end      
         
        for j=1:length(calidad)
            
            disp(['Procesando ' num2str(tiempo) ' segundos de la senal ' num2str(senales(i))...
            ' con calidad objetivo ' num2str(calidad(j)) ' y diccio ' labeldiccio]);
            ticsenal=tic;
          
            resultado(j).senal=num2str(senales(i));
            resultado(j).no_der=1;
            resultado(j).long_signal=tiempo;
            resultado(j).tipo_cal=tipo_cal;
            resultado(j).calidadobj=calidad(j);
            resultado(j).Q=Q;
            resultado(j).ultimo=0;
            
            [resultado(j).CRlow,resultado(j).CRhigh,resultado(j).prd,resultado(j).prd1,resultado(j).prdcc,resultado(j).prdseg,resultado(j).prd1seg,...
            resultado(j).CRlowseg,resultado(j).CRhighseg,resultado(j).RMSE,resultado(j).RMSEseg,resultado(j).maxidif,resultado(j).entrada,...
            resultado(j).salida,resultado(j).entradamV,resultado(j).salidamV,resultado(j).Qfinal,resultado(j).coeficientes,resultado(j).atomos]...
            =comp_MP_entrop(resultado(j).senal,resultado(j).no_der,resultado(j).long_signal,...
            resultado(j).tipo_cal,resultado(j).calidadobj,resultado(j).Q,resultado(j).ultimo,diccio);
        
            % No queremos guardar la direccion completa solo el nombre,
            % pero a la funcion hay que pasarle la  direccion completa
            resultado(j).diccio=labeldiccio;
                                
            toc(ticsenal);
        end
    
    feval('save',['senal_' resultado(j).senal '_tiempo_' num2str(tiempo) '.mat'],'resultado');
    disp(['Fin procesado senal ' resultado(j).senal]);
    toc(ticgrupo);
end

cd(oldfolder);
