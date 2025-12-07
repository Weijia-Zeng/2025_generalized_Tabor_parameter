clear; clc

lambda_E_select = [1 10 1000];

%% Horizontal axis: lambda_T

for K=1:length(lambda_E_select)
    filename1 = ['lambda_E_',num2str(lambda_E_select(K)),'_Rf_20_Rs_20'];
    load(filename1)
    Fmax1 = Fmax;
    %Fmax1(end) = NaN;
    lambda_T1 = lambda_T;
    lambda_E1 = lambda_E;

    filename2 = ['lambda_E_',num2str(lambda_E_select(K)),'_Rf_100_Rs_20'];
    load(filename2)
    Fmax2 = Fmax;
    lambda_T2 = lambda_T;
    lambda_E2 = lambda_E;

    filename3 = ['lambda_E_',num2str(lambda_E_select(K)),'_Rf_500_Rs_20'];
    load(filename3)
    Fmax3 = Fmax;
    % Fmax3(45:51) = NaN;
    lambda_T3 = lambda_T;
    lambda_E3 = lambda_E;

    f = figure(K);
    % f.Position = [1000 398 560 840];
    hold on

    plot(lambda_T1,Fmax1/pi./lambda_E1.^2,'linewidth',2,'Color','#A8D0DB')
    plot(lambda_T2,Fmax2/pi./lambda_E2.^2,'linewidth',2,'Color','#E49273')
    plot(lambda_T3,Fmax3/pi./lambda_E3.^2,'linewidth',2,'Color','#A37A74')

    xl = xline(lambda_E,'--k',['$\lambda_E=$',num2str(lambda_E_select(K))],'LineWidth',1,'Interpreter','latex','FontSize',12);
    xl.LabelVerticalAlignment = 'bottom';
    xl.LabelHorizontalAlignment = 'left';
    xl.LabelOrientation = 'horizontal';

    c = [.01 .72 .77
        .99 .49 .00
        .68 .42 .89];

    % scatter([lambda(1),lambda(8),lambda(11)],[Fmax2(1)/pi/lambda(1)^2,Fmax2(8)/pi/lambda(8)^2,Fmax2(11)/pi/lambda(11)^2],75,c,'Marker','x','LineWidth',2)
    % scatter([lambda(41),lambda(48),lambda(51)],[Fmax2(41)/pi/lambda(41)^2,Fmax2(48)/pi/lambda(48)^2,Fmax2(51)/pi/lambda(51)^2],75,c,'Marker','x','LineWidth',2)

    line([1e-3,10^(-1.5)],[2,2],'linewidth',2.5,'linestyle','--','color','k')
    line([10^(1.5),1e3],[1,1],'linewidth',2.5,'linestyle','--','color','k')

    axis on
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.MinorGridAlpha = 0.1;
    ax.GridAlpha = 0.1;
    ax.TickLength = [0.02 0.02];
    ax.XLim = [1e-3,1e3];
    ax.XScale = "log";
    ax.YLim = [0.85,2.15];
    %ax.XDir = 'reverse';
    xlabel('$\lambda_T$','interpreter','latex','fontsize',12)
    ylabel('$F/{\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)

    leg = {'$\bar{R_f}=20$','$\bar{R_f}=100$','$\bar{R_f}=500$'};

    le = legend(leg,"Interpreter","latex",'location','best','fontsize',12);
    le.EdgeColor = 'w';
    title(['$\quad \bar{R_s}=$',num2str(Rs),'$\quad \lambda_p=\infty \quad \lambda_T=\infty$'],'Interpreter','latex','fontsize',12)

end


%% horizontal axis: lambda

f = figure(K+1);
% f.Position = [1000 398 560 840];
hold on

for K=1:length(lambda_E_select)
    filename1 = ['lambda_E_',num2str(lambda_E_select(K)),'_Rf_20_Rs_20'];
    load(filename1)
    Fmax1 = Fmax;
    %Fmax1(end) = NaN;
    lambda_T1 = lambda_T;
    lambda_E1 = lambda_E;
    lambda1 = lambda;

    filename2 = ['lambda_E_',num2str(lambda_E_select(K)),'_Rf_100_Rs_20'];
    load(filename2)
    Fmax2 = Fmax;
    lambda_T2 = lambda_T;
    lambda_E2 = lambda_E;
    lambda2 = lambda;

    filename3 = ['lambda_E_',num2str(lambda_E_select(K)),'_Rf_500_Rs_20'];
    load(filename3)
    Fmax3 = Fmax;
    % Fmax3(45:51) = NaN;
    lambda_T3 = lambda_T;
    lambda_E3 = lambda_E;
    lambda3 = lambda;

    plot(lambda1,Fmax1/pi./lambda_E1.^2,'linewidth',2,'Color','#A8D0DB')
    plot(lambda2,Fmax2/pi./lambda_E2.^2,'linewidth',2,'Color','#E49273')
    plot(lambda3,Fmax3/pi./lambda_E3.^2,'linewidth',2,'Color','#A37A74')

end

line([1e-3,10^(-1.5)],[2,2],'linewidth',2.5,'linestyle','--','color','k')
line([10^(1.5),1e3],[1,1],'linewidth',2.5,'linestyle','--','color','k')

axis on
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.MinorGridAlpha = 0.1;
ax.GridAlpha = 0.1;
ax.TickLength = [0.02 0.02];
ax.XLim = [1e-3,1e3];
ax.XScale = "log";
ax.YLim = [0.85,2.15];
%ax.XDir = 'reverse';
xlabel('$\lambda$','interpreter','latex','fontsize',12)
ylabel('$F/{\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)

leg = {'$\bar{R_f}=20$','$\bar{R_f}=100$','$\bar{R_f}=500$'};

le = legend(leg,"Interpreter","latex",'location','best','fontsize',12);
le.EdgeColor = 'w';
title(['$\quad \bar{R_s}=$',num2str(Rs),'$\quad \lambda_p=\infty \quad \lambda_T=\infty$'],'Interpreter','latex','fontsize',12)
