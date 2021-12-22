
clear all
num=50;

k15=1.5; k16p=1; k16pp=2; J15=0.01; J16=0.01;

x=linspace(0,1,num);
m=linspace(0,1,num);
for i=1:num
    for j=1:num
        goldfun(i,j)=gold(k15*m(i),k16p+k16pp*x(j),J15,J16);
    end
end


surf(goldfun)
