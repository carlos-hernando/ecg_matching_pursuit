function [longlados, val, valcoef]=obtener_longlados(atomos, coeficientes, diccio)
%Obtener las longitudes de los lados de un conjunto de atomos
%que viene descrito por los indices de los ficheros en el diccionario
%y de los propios atomos en el fichero

diractual=pwd;

ficheros=dir(diccio);
cd(diccio);

h=1;

for i=1:size(atomos,1)
    if atomos(i,1)>0
        longlados(h,:)=zeros(1,2);
        feval('load',ficheros(atomos(i,1)).name);
        senalatomo=diccio.atomos(:,atomos(i,2));
        aux=find(senalatomo);
        ini=aux(1);
        fin=aux(length(aux));
        mid=find(senalatomo==max(senalatomo));
        longlados(h,1)=mid-ini+1;
        longlados(h,2)=fin-mid+1;
        h=h+1;
    end        
end


%Histograma con la cantidad de lados de diferente longitud
val=hist3(longlados, [max(longlados(:,1))-min(longlados(:,1))+1, max(longlados(:,2))-min(longlados(:,2))+1]);

%Histograma con la cantidad de energía de los atomos de las diferentes
%longitudes de lado
valcoef = zeros(size(val));
h=1;
for i=1:length(coeficientes)
    if coeficientes(i) ~= 0
        valcoef(longlados(h,1)-min(longlados(:,1))+1,longlados(h,2)-min(longlados(:,2))+1)=...
        valcoef(longlados(h,1)-min(longlados(:,1))+1,longlados(h,2)-min(longlados(:,2))+1)+coeficientes(i)^2;
        h=h+1;
    end
end

dibujar_histogramas(longlados,val,valcoef);

cd(diractual);


