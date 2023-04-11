function stateplt(time,yo,xt,x,cvar)
figure('Units','inches','Position',[5 1 9.5 6.5]);clf
subplot(2,1,1)
plot(time,xt,'r');hold on
plot(time,x,'b');
plot(time,yo,'k*','MarkerSize',3);
h=legend('Tru','Ana');
set(gca,'Position',[0.075 0.6097 0.8 0.28])
set(h,'Position',[0.88 0.8 0.1 0.1])
legend('boxoff')
xlabel('Time')
ylabel(cvar,'Fontsize',12)
title('Truth vs. Analysis','Fontsize',12,'Fontweight','bold')

subplot(2,1,2)
plot(time,x-xt,'b','Linewidth',1.5);hold on
plot([time(1) time(end)],[0 0],'-k')
set(gca,'Ylim',[-2 2],'Xlim',[time(1) time(end)])
set(gca,'Position',[0.075 0.135 0.8 0.28])
title('Analysis-Truth','Fontsize',12,'Fontweight','bold')
xlabel('Time')
ylabel('Error')