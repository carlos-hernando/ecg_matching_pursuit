function[prob,entrop,simbolos]=calc_entrop(secuencia)

%Calcula la entropía de una secuencia de símbolos

cuenta(1)=1;
simbolos(1)=secuencia(1);

for i=2:length(secuencia)
    flag=0;
    for j=1:length(simbolos)
        if secuencia(i)==simbolos(j)
            cuenta(j)=cuenta(j)+1;
            flag=1;
            break;
        end
    end
    if flag==0
        cuenta(j+1)=1;
        simbolos(j+1)=secuencia(i);
    end
end
            
prob=cuenta/length(secuencia);

entrop=sum(-prob.*log2(prob));