%
set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);


%specify parameters
vi=0.1
ve=0.5;
vt=vi+ve;
st=1;
si=[0.01:.01:10];
X=(si.^2*(ve^2-vi^2)+2*si*st*vi*vt-st^2*vt^2)./si;
figure(1)
plot(X,si,'linewidth',2)

axis([0 10 0 10])

