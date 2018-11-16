%
%	Cálculo de la S/N según fórmula de Nguyen
%	(Wavelet and filters banks, pág.385)
%
%	Parámetros de entrada:
%		xorg=señal original
%		xrec=señal reconstruida
%
%	Parámetros de salida:
%		prd= raíz de la diferencia cuadrática media 
%

function prd=prdmedio(xorg,xrec)
if nargin==2
	if length(xorg)==length(xrec)
      difer=xorg-xrec;
      xorgmedio=mean(xorg);
		prd=sqrt(sum((difer).^2)/sum((xorg-xorgmedio).^2))*100;
	else
		error('la longitud de la señales original y reconstruida debe ser la misma')
	end
else
   error('el número de argumentos debe ser dos')
end