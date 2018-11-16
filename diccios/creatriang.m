function triangulo=creatriang(iz,de,altura)
% Funcion que crea un triangulo segundo las longitudes de los lados y la
% altura indicada

TL=1;
TC=iz+1;
TR=iz+de+1;

triangulo=zeros(1,TR);
triangulo(1:TC-1) = altura*([1:TC-1]-TL)/(TC-TL);
triangulo(TC+1:TR)= altura*(TR-[TC+1:TR])/(TR-TC);
triangulo(TC)=altura;