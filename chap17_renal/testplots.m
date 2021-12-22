
tau=4; gamma=4;

test=@(x) sin(x*tau+x) - sin(x*tau);

fplot(@(x) test(x), [0,15])