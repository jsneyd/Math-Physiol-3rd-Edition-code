%here I am trting out som stocahstic simulations.
clear all

dt = 0.1;

k(1,1) = 0;
k(1,2) = 1;
k(1,3) = 0.6;


k(2,1) = 0.4;
k(2,2) = 0;
k(2,3) = .1;

k=dt*k;
%Form the test matrix

K = sum(k,2);

k(1,1) = 1-K(1);
k(2,2) = 1-K(2);

Tst = cumsum(k,2);


s = 1;
S(1) = s;
j = 1;
while(s<3)
rtst = rand(1);
    news = min(find(Tst(s,:)>=rtst));
    j=j+1;
    s = news;
    S(j) = s;
end;

    
    