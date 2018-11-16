function y_i=successive_convolutions(L,ith)

x=ones(1,L);
y_i=x;
for k=1:ith
    y_i=(conv(y_i,x));
end
y_i=y_i/max(y_i);