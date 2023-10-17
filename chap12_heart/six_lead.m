%heart vectors

I1=[1,0];
I2=[1/2,-sqrt(3)/2];
I3=[-1/2,-sqrt(3)/2];
Iavf=[0,-1];
Iavr=[-sqrt(3)/2,1/2];
Iavl= [sqrt(3)/2,1/2];

A=[I1;I2;I3;Iavr;Iavl;Iavf]



R=[7;9;0;-8;3;7];

B=A'*A;
Rb = A'*R;
x=B\Rb