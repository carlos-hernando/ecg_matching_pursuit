function [param]=indicetoparam(indice)

%Esta funcion proprociona los parametros de longitud de
%señal cuadrada y número de convoluciones que le corresponde a
%un determinado medio atomo descrito por su indice unívoco 
%
%   Parámetros de entrada:
%   indice: indice univoco que define al atomo
%
%   Parámetros de salida
%   param: vector 2x1 con los datos [1,1]=longitud señal cuadrada,
%   [2,1]=numero de convoluciones
%              ______________________________________________
%             |                                              |
%             | ULTIMA ACTUALIZACION: 18 de octubre 2014     |
%             | AUTOR: Carlos Hernando Ramiro                |
%             |______________________________________________|
%

iteracion=1;
stop=0;
D=2;

if (isempty(indice))
    disp('Error indice no valido');
    param=[];
else
    %Si indice=1 entonces es el caso especial de 0 convoluciones, que puee ser
    %parte de un atomodelta o parte ausente de medio atomo
    if indice == 1
        param=[0;0];
    else
        while(stop~=1)
            LTe=2*D;
            LTo=2*D-1;
            stage=1;
            while ((fix((LTe+stage)/(stage+1))>=2 || fix((LTo+stage)/(stage+1))>=2) && stop~=1)
                if rem((LTe+stage),(stage+1))==0
                    iteracion=iteracion+1;
                    if iteracion == indice
                        param=[(LTe+stage)/(stage+1) ; stage];
                        stop=1;
                    end
                end
                if rem((LTo+stage),(stage+1))==0
                    iteracion=iteracion+1;
                    if iteracion == indice
                        param=[(LTo+stage)/(stage+1) ; stage];
                        stop=1;
                    end
                end
                stage=stage+1;
            end
            D=D+1;
        end
    end
end        
        