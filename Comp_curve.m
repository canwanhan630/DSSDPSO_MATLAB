% clear all;
% close all;
% clc;

T = 1:500;

figure;
hold on
plot(T,DSSDPSO,'b');
plot(T,NPSOWM,'-*r');
plot(T,RDPSO,'-ok');
set(gca,'FontSize',15);
xlabel('��������','fontsize',15)
ylabel('��ֵ�԰��ƽ(dB)','fontsize',15)
legend('DSSDPSO','NPSOWM','RDPSO','fontsize',13)
hold off