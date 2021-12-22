a = 3/4;

x = [0:.01:100];
V  = -x.*(erfc(sqrt(x)*a)-erfc(sqrt(x)*(1-a)));
plot(x,V,'linewidth',2)
xlabel('x','fontsize',16)
ylabel('v(x)','fontsize',16)

out=[x' V'];

save('flasher.out','-ascii','out')