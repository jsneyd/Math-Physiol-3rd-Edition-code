
load clotting

scale=50;

figure(1)

plot(u1(20/delt,:),'blue')
hold on
plot(u3(20/delt,:)/scale,'--blue')
plot(u1(100/delt,:),'red')
plot(u3(100/delt,:)/scale,'--red')
plot(u1(150/delt,:),'green')
plot(u3(150/delt,:)/scale,'--green')
plot(u1(200/delt,:),'black')
plot(u3(200/delt,:)/scale,'--black')
hold off

dum=[u1(20/delt,:)' u1(100/delt,:)' u1(150/delt,:)' u1(200/delt,:)' u3(20/delt,:)' u3(100/delt,:)' u3(150/delt,:)' u3(200/delt,:)'];
save('clotting.dat','-ascii','dum')


% for i=1:35
%     filteru1(i,:)=u1(100*i,:);
%     filteru2(i,:)=u2(100*i,:);
%     filteru3(i,:)=u3(100*i,:);
% end
% 
% figure(1)
% surface(filteru1)
% 
% figure(2)
% surface(filteru2)
% 
% figure(3)
% surface(filteru3)