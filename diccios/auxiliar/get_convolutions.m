function L_stage=get_convolutions(D)

LTe=2*D;
LTo=2*D-1;
stage=1;
L=[];
n_conv=[];
while fix((LTe+stage)/(stage+1))>=2||fix((LTo+stage)/(stage+1))>=2
    if rem((LTe+stage),(stage+1))==0
        L=[L (LTe+stage)/(stage+1)];
        n_conv=[n_conv stage];
    end
    if rem((LTo+stage),(stage+1))==0
        L=[L (LTo+stage)/(stage+1)];
        n_conv=[n_conv stage];
    end
    stage=stage+1;
end
L_stage=[L; n_conv];