function [indice]=paramtoindice(param)

%Esta funcion proprociona el indice unívoco que le corresponde a
%un determinado medio atomo definido por sus parametros de longitud de
%señal cuadrada y número de convoluciones
%
%   Parámetros de entrada:
%   param: vector 2x1 con los datos [1,1]=longitud señal cuadrada,
%   [2,1]=numero de convoluciones
%
%   Parámetros de salida
%   indice: indice univoco que define al atomo
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 18 de octubre 2014     |
%             | AUTOR: Carlos Hernando Ramiro                |
%             |______________________________________________|
%

indice=[];
stop=0;

%Usamos la misma nomenclatura que en el programa get_convolutions.m de
%Manuel Blanco
Dtarget=param(1);
stagetarget=param(2);
tamano=ceil((((stagetarget+1)*Dtarget)-stagetarget)/2);
D=2;

indice=1;
%Si se obrepasa el tamano sin stop los param no son validos o es delta
while(stop~=1 && D<=tamano)
    LTe=2*D;
    LTo=2*D-1;
    stage=1;
    while ((fix((LTe+stage)/(stage+1))>=2 || fix((LTo+stage)/(stage+1))>=2) && stop~=1)
        if rem((LTe+stage),(stage+1))==0
            indice=indice+1;
            if ((LTe+stage)/(stage+1))==Dtarget && stage==stagetarget
                stop=1;
            end
        end
        if rem((LTo+stage),(stage+1))==0
            indice=indice+1;
            if ((LTo+stage)/(stage+1))==Dtarget && stage==stagetarget
                stop=1;
            end
        end
        stage=stage+1;
    end
    D=D+1;
end

%Si no es delta no tiene indice
if indice==1 && stagetarget~=0 && Dtarget~=0
    disp('Error parametros no validos');
    indice=[];
end
