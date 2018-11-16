%
% Max error es una funci�n que calcula error m�ximo entre dos se�ales
% 	seg�n f�rmula  propuesta en Nguyen (Wavelet and filters banks, p�g.385)
%
%	Par�metros de entrada:
%		xor=se�al original
%		xr=se�al reconstruida
%
%	Par�metros de salida:
%		valor= maximo error entre las dos se�ales 
%

function valor=maxerror(xorg,xrec)
if nargin==2
	if length(xorg)==length(xrec)
		difer=xorg-xrec;
		valor=max(abs(difer));
	else
		error('la longitud de la se�ales original y reconstruida debe ser la misma')
	end
else
	error('el n�mero de argumentos debe ser dos')
end