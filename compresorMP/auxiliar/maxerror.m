%
% Max error es una función que calcula error máximo entre dos señales
% 	según fórmula  propuesta en Nguyen (Wavelet and filters banks, pág.385)
%
%	Parámetros de entrada:
%		xor=señal original
%		xr=señal reconstruida
%
%	Parámetros de salida:
%		valor= maximo error entre las dos señales 
%

function valor=maxerror(xorg,xrec)
if nargin==2
	if length(xorg)==length(xrec)
		difer=xorg-xrec;
		valor=max(abs(difer));
	else
		error('la longitud de la señales original y reconstruida debe ser la misma')
	end
else
	error('el número de argumentos debe ser dos')
end