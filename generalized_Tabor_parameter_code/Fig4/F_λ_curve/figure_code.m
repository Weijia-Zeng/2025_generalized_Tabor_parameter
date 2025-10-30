clear; clc

load('Rf_20_Rs_20.mat')
Fmax1 = Fmax;
Fmax1(end) = NaN;
lambda_E1 = lambda_E;

load('Rf_100_Rs_20.mat')
Fmax2 = Fmax;
lambda_E2 = lambda_E;

load('Rf_500_RS_20.mat')
Fmax3 = Fmax;
Fmax3(45:51) = NaN;
lambda_E3 = lambda_E;

leg = {'$\bar{R_f}=20$','$\bar{R_f}=100$','$\bar{R_f}=500$'};


f = figure(2);
% f.Position = [1000 398 560 840];
hold on

plot(lambda_E1,Fmax1/pi./lambda_E1.^2,'linewidth',2,'Color','#A8D0DB')
plot(lambda_E2,Fmax2/pi./lambda_E2.^2,'linewidth',2,'Color','#E49273')
plot(lambda_E3,Fmax3/pi./lambda_E3.^2,'linewidth',2,'Color','#A37A74')

line([1e-3,10^(-1.5)],[2,2],'linewidth',2.5,'linestyle','--','color','k')
line([10^(0.5),1e2],[1,1],'linewidth',2.5,'linestyle','--','color','k')

c = parula(10);

scatter([lambda(1),lambda(11),lambda(21)],[Fmax2(1)/pi/lambda(1)^2,Fmax2(11)/pi/lambda(11)^2,Fmax2(21)/pi/lambda(21)^2],75,c(3:5,:),'Marker','x','LineWidth',2)
scatter([lambda(31),lambda(41),lambda(51)],[Fmax2(31)/pi/lambda(31)^2,Fmax2(41)/pi/lambda(41)^2,Fmax2(51)/pi/lambda(51)^2],75,c(6:8,:),'Marker','x','LineWidth',2)

axis on
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.MinorGridAlpha = 0.1;
ax.GridAlpha = 0.1;
ax.TickLength = [0.02 0.02];
ax.XLim = [1e-3,1e2];
ax.XScale = "log";
ax.YLim = [0.85,2.15];
%ax.XDir = 'reverse';
xlabel('$\lambda_E$','interpreter','latex','fontsize',12)
ylabel('$F/{\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)
le = legend(leg,"Interpreter","latex",'location','best','fontsize',12);
le.EdgeColor = 'w';
title(['$\quad \bar{R_s}=$',num2str(Rs),'$\quad \lambda_p=\infty \quad \lambda_T=\infty$'],'Interpreter','latex','fontsize',12)