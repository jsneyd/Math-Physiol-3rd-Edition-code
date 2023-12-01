%stability for Goldwin model
p=[8.1:.01:20];

y=(8./(p-8)).^(1./p);

b = 1./((1+y.^p).*y);


figure(1)
plot(b,p)
 
