%
%	C�lculo de la S/N seg�n f�rmula de Nguyen
%	(Wavelet and filters banks, p�g.385)
%
%	Par�metros de entrada:
%		xorg=se�al original
%		xrec=se�al reconstruida
%
%	Par�metros de salida:
%		prd= ra�z de la diferencia cuadr�tica media 
%

function prd=prdmedio(xorg,xrec)
if nargin==2
	if length(xorg)==length(xrec)
      difer=xorg-xrec;
      xorgmedio=mean(xorg);
		prd=sqrt(sum((difer).^2)/sum((xorg-xorgmedio).^2))*100;
	else
		error('la longitud de la se�ales original y reconstruida debe ser la misma')
	end
else
   error('el n�mero de argumentos debe ser dos')
end