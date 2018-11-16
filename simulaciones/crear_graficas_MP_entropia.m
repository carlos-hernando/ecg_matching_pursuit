function [prdmed, CRlowmed, CRhighmed, RMSEmed, NSAmed, MAXmed]=crear_graficas_MP_entropia(carpeta)

diractual = pwd;

ficheros=dir(carpeta);
cd(carpeta);

i=3;


while i<=length(ficheros)
    
       
        [prd, CRlow, CRhigh, RMSE, NSA, MAX]=tomar_datos([carpeta '\' ficheros(i).name]);
        
        
        prdtod(i-2,:, :) = prd;
        CRlowtod(i-2,:, :) = CRlow;
        CRhightod(i-2,:, :) = CRhigh;
        RMSEtod(i-2,:, :) = RMSE;
        NSAtod(i-2,:, :) = NSA;
        MAXtod(i-2,:, :) = MAX;
        
        ref = strfind(ficheros(i).name,'dic');
        ref2 = strfind(ficheros(i).name,'tip');
        s{i-2}= regexprep(ficheros(i).name(ref:ref2-1), '_', ' ');     
    
        i = i+1;
end

for i=1:size(prdtod,1)
    
    prdmed(i,:) = mean(prdtod(i,:,:));
    CRlowmed(i,:) = mean(CRlowtod(i,:,:));
    CRhighmed(i,:) = mean(CRhightod(i,:,:));
    RMSEmed(i,:) = mean(RMSEtod(i,:,:));
    NSAmed(i,:) = mean(NSAtod(i,:,:));
    MAXmed(i,:) = mean(MAXtod(i,:,:));
    
end

h = zeros(1, size(prdmed,1));

colors = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 1; 1 0 0; 1 1 0; 0.5 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5; 1 0.5 0.5; 0 0 0.5; 0 0.5 0; 0.5 0 0; 0.5 0.5 0];

%CRlow results:

figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('CRlow')
hold on;

for k=1:size(prdmed,1)
    h(k) = plot(prdmed(k,:), CRlowmed(k,:), '-o', 'Color', colors(mod(k-1,length(colors))+1,:), 'LineWidth',2,...
    'MarkerEdgeColor',colors(mod(k-1,length(colors))+1,:),...
    'MarkerFaceColor', colors(mod(k-1,length(colors))+1,:),'MarkerSize',6);
end

legend(h, s, 'Location', 'NorthWest');


%CRhigh results:

figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('CRhigh')
hold on;

for k=1:size(prdmed,1)
    h(k) = plot(prdmed(k,:), CRhighmed(k,:), '-o', 'Color', colors(mod(k-1,length(colors))+1,:), 'LineWidth',2,...
    'MarkerEdgeColor',colors(mod(k-1,length(colors))+1,:),...
    'MarkerFaceColor', colors(mod(k-1,length(colors))+1,:),'MarkerSize',6);
end

legend(h, s, 'Location', 'NorthWest');


%NSA results

figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('NSA')
hold on;

for k=1:size(prdmed,1)
    h(k) = plot(prdmed(k,:), NSAmed(k,:), '-s', 'Color', colors(mod(k-1,length(colors))+1,:), 'LineWidth',2,...
    'MarkerEdgeColor',colors(mod(k-1,length(colors))+1,:),...
    'MarkerFaceColor', colors(mod(k-1,length(colors))+1,:),'MarkerSize',6);
end

legend(h, s, 'Location', 'NorthWest');


%MAX results
figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('MAX (\muV)')
hold on;

for k=1:size(prdmed,1)
    h(k) = plot(prdmed(k,:), MAXmed(k,:).*1000, '-^', 'Color', colors(mod(k-1,length(colors))+1,:), 'LineWidth',2,...
    'MarkerEdgeColor',colors(mod(k-1,length(colors))+1,:),...
    'MarkerFaceColor', colors(mod(k-1,length(colors))+1,:),'MarkerSize',6);
end

legend(h, s, 'Location', 'NorthWest');

cd(diractual);
