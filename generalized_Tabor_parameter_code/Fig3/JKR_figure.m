clear; clc

%% F-delta curve
figure(1)

figure_width = 8;
figure_height = 8*0.75;
figure_number_fontsize = 10;
line_width = 1.5;

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

Index = zeros(1,3);

hold on

c = parula(4);
leg = cell(1,3);

for i = 1:3
    switch i
        case 1
            load('JKR_limit_P_5e-05_T_0.01_Gamma_0.0001_alpha_0.05.mat')
        case 2
            load('JKR_limit_P_5e-05_T_0.01_Gamma_0.01_alpha_0.05.mat')
        case 3
            load('JKR_limit_P_5e-05_T_0.01_Gamma_1_alpha_0.05.mat')
    end

    [Fmax,Index(i)] = max(RF1);

    plot((Rd(Index(i):end)-delta0),RF1(Index(i):end),'Color',c(i,:),'LineStyle','-','LineWidth',line_width)

    % if i == 1
    %     leg{i} = ['$P=$',num2str(P,'%.5f')];
    % else
         leg{i} = ['$\Gamma=$',num2str(Gamma)];
    % end

end

for i = 1:3
    switch i
        case 1
            load('JKR_limit_P_5e-05_T_0.01_Gamma_0.0001_alpha_0.05.mat')
        case 2
            load('JKR_limit_P_5e-05_T_0.01_Gamma_0.01_alpha_0.05.mat')
        case 3
            load('JKR_limit_P_5e-05_T_0.01_Gamma_1_alpha_0.05.mat')
    end

    plot((Rd(1:Index(i))-delta0),RF1(1:Index(i)),'Color',c(i,:),'LineStyle','--','LineWidth',line_width)
end

% y_line = yline(0,'-k','alpha',0.8);

axis on
box on
ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.MinorGridAlpha = 0.1;
ax.GridAlpha = 0.1;
ax.TickLength = [0.02 0.02];
set(gca,'FontSize',figure_number_fontsize)
ax.YLim = [-0.1,1.1];
legend(leg,'Interpreter','latex','location','best','EdgeColor','none','Color','White','FontSize',8)
xlabel('$\bar{\Delta}=\delta/\delta_*$','interpreter','latex')
ylabel('$\mathcal{F}=F/(\pi\gamma {R_s})$','interpreter','latex')

% zp = BaseZoom([8.5,0,5,0.4],[-4,-0.1,3,0.2]);
% zp.run
% zp.subAxes.Color = 'white';
% zp.subAxes.YTick = [0];

hold off


%% F-A curve
figure(2)

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

c = parula(4);

for i = 1:3
    switch i
        case 1
            load('JKR_limit_P_5e-05_T_0.01_Gamma_0.0001_alpha_0.05.mat')
        case 2
            load('JKR_limit_P_5e-05_T_0.01_Gamma_0.01_alpha_0.05.mat')
        case 3
            load('JKR_limit_P_5e-05_T_0.01_Gamma_1_alpha_0.05.mat')
    end

    plot(RA(Index(i):end),RF1(Index(i):end),'Color',c(i,:),'LineStyle','-','LineWidth',line_width)
end

for i = 1:3
    switch i
        case 1
            load('JKR_limit_P_5e-05_T_0.01_Gamma_0.0001_alpha_0.05.mat')
        case 2
            load('JKR_limit_P_5e-05_T_0.01_Gamma_0.01_alpha_0.05.mat')
        case 3
            load('JKR_limit_P_5e-05_T_0.01_Gamma_1_alpha_0.05.mat')
    end

    plot(RA(1:Index(i)),RF1(1:Index(i)),'Color',c(i,:),'LineStyle','--','LineWidth',line_width)
end

axis on
box on
ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.MinorGridAlpha = 0.1;
ax.GridAlpha = 0.1;
ax.TickLength = [0.02 0.02];
set(gca,'FontSize',figure_number_fontsize)
ax.YLim = [-0.1,1.1];
legend(leg,'Interpreter','latex','location','best','EdgeColor','none','Color','White','FontSize',8)
xlabel('$A=a/a_*$','interpreter','latex')
ylabel('$\mathcal{F}=F/(\pi\gamma {R_s})$','interpreter','latex')

hold off




