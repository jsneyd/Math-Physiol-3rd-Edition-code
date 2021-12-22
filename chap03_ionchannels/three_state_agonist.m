
% initialize the states
clear
% this simulates the agonist binding channel model
%set the rate constants
r(1) = 2; % alpha  (close rate)
r(2) = 10; %gamma  (binding rate)
r(3) = 1;  %beta  (open rate)
r(4) = 0.5;  %delta (unbinding rate)
N = 4;
T = [];
t = [];
xx = [];
X = [];

x = 2;
T(1) = 0;
t(1) = 0;
X(1) = x;
xx(1) = x;
j = 2;
jj = 2;
jo = 1;
jc = 1;
oldj = 1;
ot = [];
ct = [];


R(1,1) = r(1);  %1->2 (close channel)
R(2,1) = 0;  %null transtion
R(1,2) = r(3); % 2->1 (open channel)
R(2,2) = r(4);  %2->3  (unbinding)
R(1,3) = r(2); %3->2(binding)
R(2,3) = 0;

%the transitions:
C(1,1) = 2;
C(1,2) = 1;
C(2,2) =3;
C(1,3) = 2;

h = sum(R);

oldstate = x;
T(1) = 0;
X(1) = x;
j = 2;
for j = 2:5000
rn = rand(1);
    H = h(x);
    dt = -log(rn)/H;
    T(j) =  dt+T(j-1);
    
tm = cumsum(R(:,x))/H;
rm = rand(1);
    rk = min(find(rm<=tm))  ;
 xx(jj) = x;
x = C(rk,x);
X(j) = x;

t(jj) = T(j);
t(jj+1) = T(j);
xx(jj+1) = x;
jj = jj+2;

if ((oldstate==2)*(x==1))  %channel opened
ct(jc) = T(j) - T(oldj);
oldstate = 1;
jc = jc + 1;
oldj = j;
end

if(oldstate==1)  %channel closed
    ot(jo) = dt;
    jo = jo+1;
    oldj = j;
    oldstate = 2;
end

    
    
end
figure(1)

plot(t,xx)

figure(2)
hist(log(ct),jc/15)
figure(3)
hist(log(ot),jo/15)



 