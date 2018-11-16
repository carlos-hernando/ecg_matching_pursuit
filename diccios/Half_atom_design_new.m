function[Half_atom, L_stage]=Half_atom_design_new(D)
% Funcion igual que Half_atom_design de Manuel Blanco pero que devuelve 
% tambien los parametros de convolucion

L_stage=get_convolutions(D);
Half_atom=zeros(size(L_stage,2),D+1); %we add +1 to include a 0
for k=1:size(L_stage,2)
    y=successive_convolutions(L_stage(1,k),L_stage(2,k));
    Half_atom(k,:)=[0 y(1:D)];
end