function[]=dibujar_histogramas(longlados,val,valcoef)

%Histograma con la cantidad de lados de diferente longitud
x=min(longlados(:,1)):1:max(longlados(:,1));
y=min(longlados(:,2)):1:max(longlados(:,2));
figure;
imagesc(x,y,val');
colorbar;
set(gca,'XTick',[min(x):max(x)]);
set(gca,'YTick',[min(y):max(y)]);
xlabel('Left-half length');
ylabel('Right-half length');
title('Number of atoms');

%Histograma con la cantidad de energía de los atomos de las diferentes
%longitudes de lado
figure;
imagesc(x,y,valcoef');
colorbar;
set(gca,'XTick',[min(x):max(x)]);
set(gca,'YTick',[min(y):max(y)]);
xlabel('Left-half length');
ylabel('Right-half length');
title('Energy of the atoms');