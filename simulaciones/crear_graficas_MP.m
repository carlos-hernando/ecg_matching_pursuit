function [prdmed, CR1med, CR2med, CR3med, RMSEmed, NSAmed, MAXmed]=crear_graficas_MP(carpeta)

diractual = pwd;

ficheros=dir(carpeta);
cd(carpeta);

i=3;


while i<=length(ficheros)

    [prd, CR, CR1, CR2, CR3, RMSE, NSA, MAX]=tomar_datos([carpeta '\' ficheros(i).name]);

    prdtod(i-2,:, :) = prd;
    CR1tod(i-2,:, :) = CR1;
    CR2tod(i-2,:, :) = CR2;
    CR3tod(i-2,:, :) = CR3;
    CRtod(i-2,:, :) = CR;
    RMSEtod(i-2,:, :) = RMSE;
    NSAtod(i-2,:, :) = NSA;
    MAXtod(i-2,:, :) = MAX;

    ref = strfind(ficheros(i).name,'dic');
    if isempty(ref)    
        ref=1;
        ref2 = strfind(ficheros(i).name,'bf');
        ref3 = strfind(ficheros(i).name,'seg');
        ref4 = strfind(ficheros(i).name,'_tip');
        s{i-2}= regexprep([ficheros(i).name(ref:ref2-2) ficheros(i).name(ref3+3:ref4-1)], '_', ' ');       
    else
        ref2 = strfind(ficheros(i).name,'tip');
        %Quitamos el segundo 1 ya que se emplea siempre el atomo delta
        auxname = ficheros(i).name(ref+4:ref2-2);
        auxname = [auxname(1:2) auxname(5:length(auxname))];
        s{i-2}= regexprep(auxname, '_', ' ');
    end

    i = i + 1;
    
end

for i=1:size(prdtod,1)
    
    if size(prdtod,2) > 1
        
        prdmed(i,:) = mean(prdtod(i,:,:));
        CR(i,:) = mean(CRtod(i,:,:));
        CR1med(i,:) = mean(CR1tod(i,:,:));
        CR2med(i,:) = mean(CR2tod(i,:,:));
        CR3med(i,:) = mean(CR3tod(i,:,:));
        RMSEmed(i,:) = mean(RMSEtod(i,:,:));
        NSAmed(i,:) = mean(NSAtod(i,:,:));
        MAXmed(i,:) = mean(MAXtod(i,:,:));

    else

        prdmed(i,:) = prdtod(i,:,:);
        CR(i,:) = CRtod(i,:,:);
        CR1med(i,:) = CR1tod(i,:,:);
        CR2med(i,:) = CR2tod(i,:,:);
        CR3med(i,:) = CR3tod(i,:,:);
        RMSEmed(i,:) = RMSEtod(i,:,:);
        NSAmed(i,:) = NSAtod(i,:,:);
        MAXmed(i,:) = MAXtod(i,:,:);
        
    end

end

colors = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 1; 1 0 0; 1 1 0; 0.5 0.5 0.5; 0.5 0.5 1; 0.5 1 0.5; 1 0.5 0.5; 0 0 0.5; 0 0.5 0; 0.5 0 0; 0.5 0.5 0];


%CR all methods results together:
%{
figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
%title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('CR')
hold on;
h=[];
m=1;
sfin=[];

for k=1:size(prdmed,1)
    
        h(m) = plot(prdmed(k,:), CR1med(k,:), '-o', 'Color', colors(mod(m-1,length(colors))+1,:), 'LineWidth',2,...
        'MarkerEdgeColor',colors(mod(m-1,length(colors))+1,:),...
        'MarkerFaceColor', colors(mod(m-1,length(colors))+1,:),'MarkerSize',6);
 
        sfin{m}=[s{k} ' direct'];

        m=m+1;
    
        h(m) = plot(prdmed(k,:), CR2med(k,:), '-o', 'Color', colors(mod(m-1,length(colors))+1,:), 'LineWidth',2,...
        'MarkerEdgeColor',colors(mod(m-1,length(colors))+1,:),...
        'MarkerFaceColor', colors(mod(m-1,length(colors))+1,:),'MarkerSize',6);
    
        sfin{m}=[s{k} ' Arith'];

        m=m+1;
    
        h(m) = plot(prdmed(k,:), CR3med(k,:), '-o', 'Color', colors(mod(m-1,length(colors))+1,:), 'LineWidth',2,...
        'MarkerEdgeColor',colors(mod(m-1,length(colors))+1,:),...
        'MarkerFaceColor', colors(mod(m-1,length(colors))+1,:),'MarkerSize',6);
    
        sfin{m}=[s{k} ' Huff'];

        m=m+1;
        
end

legend(h, sfin, 'Location', 'NorthWest');
%}


%CR methods individually:
%{
%Method 1
figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('CR')
hold on;
h=[];
m=1;
sfin=[];

for k=1:size(prdmed,1)
        h(m) = plot(prdmed(k,:), CR1med(k,:), '-o', 'Color', colors(mod(k-1,length(colors))+1,:), 'LineWidth',2,...
        'MarkerEdgeColor',colors(mod(k-1,length(colors))+1,:),...
        'MarkerFaceColor', colors(mod(k-1,length(colors))+1,:),'MarkerSize',6);
    
        sfin{m}=[s{k} ' direct'];
    
    m=m+1;
end

legend(h, sfin, 'Location', 'NorthWest');


%Method 2
figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('CR')
hold on;
h=[];
m=1;
sfin=[];

for k=1:size(prdmed,1)
        h(m) = plot(prdmed(k,:), CR2med(k,:), '-o', 'Color', colors(mod(k-1,length(colors))+1,:), 'LineWidth',2,...
        'MarkerEdgeColor',colors(mod(k-1,length(colors))+1,:),...
        'MarkerFaceColor', colors(mod(k-1,length(colors))+1,:),'MarkerSize',6);
    
        sfin{m}=[s{k} ' Arith'];
    
    m=m+1;
end

legend(h, sfin, 'Location', 'NorthWest');
%}

%Method 3

figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD (%)')
ylabel('CR')
%ylabel('MAX (\muV)')
hold on;
h=[];
m=1;
sfin=[];

set1=[3 1 7 5 11 9 19 17];

for i=1:length(set1)
       
    k = set1(i);
    
    h(m) = plot(prdmed(k,:), CR3med(k,:), '-o', 'Color', colors(mod(m-1,length(colors))+1,:), 'LineWidth',2,...
    'MarkerEdgeColor',colors(mod(m-1,length(colors))+1,:),...
    'MarkerFaceColor', colors(mod(m-1,length(colors))+1,:),'MarkerSize',6);

    %h(m) = plot(prdmed(k,:), MAXmed(k,:).*1000, '-^', 'Color', colors(mod(m-1,length(colors))+1,:), 'LineWidth',2,...
    %'MarkerEdgeColor',colors(mod(m-1,length(colors))+1,:),...
    %'MarkerFaceColor', colors(mod(m-1,length(colors))+1,:),'MarkerSize',6);

    sfin{m}=[s{k} ' Huff'];
    %sfin{m}=['D' num2str(i)];
    
    m=m+1;
        
end

legend(h, sfin, 'Location', 'NorthWest');
figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD (%)')
ylabel('CR')
%ylabel('MAX (\muV)')
hold on;
h=[];
m=1;
sfin=[];
set2=set1 + 1;

for i=1:length(set2)
       
    k = set2(i);
    
    h(m) = plot(prdmed(k,:), CR3med(k,:), '-o', 'Color', colors(mod(m-1,length(colors))+1,:), 'LineWidth',2,...
    'MarkerEdgeColor',colors(mod(m-1,length(colors))+1,:),...
    'MarkerFaceColor', colors(mod(m-1,length(colors))+1,:),'MarkerSize',6);

    %h(m) = plot(prdmed(k,:), MAXmed(k,:).*1000, '-^', 'Color', colors(mod(m-1,length(colors))+1,:), 'LineWidth',2,...
    %'MarkerEdgeColor',colors(mod(m-1,length(colors))+1,:),...
    %'MarkerFaceColor', colors(mod(m-1,length(colors))+1,:),'MarkerSize',6);

    sfin{m}=[s{k} ' Huff'];
    %sfin{m}=['D' num2str(i)];

    m=m+1;
        
end

legend(h, sfin, 'Location', 'NorthWest');
%}


%NSA results
%{
figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('NSA')
hold on;
h=[];
m=1;
sfin=[];

for k=1:size(prdmed,1)
    if NSAmed(k,1) > 0
        h(m) = plot(prdmed(k,:), NSAmed(k,:), '-s', 'Color', colors(mod(k-1,length(colors))+1,:), 'LineWidth',1,...
        'MarkerEdgeColor',colors(mod(k-1,length(colors))+1,:),...
        'MarkerFaceColor', colors(mod(k-1,length(colors))+1,:),'MarkerSize',6);
        
        sfin{m}=s{k};
        
        m=m+1;
    end
    
    
end

legend(h, sfin, 'Location', 'NorthWest');
%}
%{
%MAX results
figure;
ref = strfind(carpeta,'senal');
ref2 =length(carpeta);
title(regexprep(carpeta(ref:ref2), '_', ' '));
xlabel('PRD')
ylabel('MAX (\muV)')
hold on;
h=[];
m=1;

for k=1:size(prdmed,1)
    h(m) = plot(prdmed(k,:), MAXmed(k,:).*1000, '-^', 'Color', colors(mod(k-1,length(colors))+1,:), 'LineWidth',2,...
    'MarkerEdgeColor',colors(mod(k-1,length(colors))+1,:),...
    'MarkerFaceColor', colors(mod(k-1,length(colors))+1,:),'MarkerSize',6);

    m=m+1;
end

legend(h, s, 'Location', 'NorthWest');
%}
cd(diractual);
