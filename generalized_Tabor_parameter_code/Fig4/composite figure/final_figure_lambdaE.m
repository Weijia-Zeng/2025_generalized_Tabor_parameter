clear; clc

figure(1)

figure_width = 8;
figure_height = 12;
figure_number_fontsize = 10;
line_width = 1.5;

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

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

% plot(lambda_E1,Fmax1/pi./lambda_E1.^2,'linewidth',2,'Color','#A8D0DB')
% plot(lambda_E2,Fmax2/pi./lambda_E2.^2,'linewidth',2,'Color','#E49273')
plot(lambda_E2,Fmax2/pi./lambda_E2.^2,'linewidth',line_width,'Color','k')
% plot(lambda_E3,Fmax3/pi./lambda_E3.^2,'linewidth',2,'Color','#A37A74')

line([1e-3,10^(-1)],[2,2],'linewidth',2.5,'linestyle','--','color','k')
line([10^(0),1e2],[1,1],'linewidth',2.5,'linestyle','--','color','k')

c = parula(10);

scatter([lambda(1),lambda(11),lambda(21)],[Fmax2(1)/pi/lambda(1)^2,Fmax2(11)/pi/lambda(11)^2,Fmax2(21)/pi/lambda(21)^2],125,c(3:5,:),'Marker','x','LineWidth',2)
scatter([lambda(31),lambda(41),lambda(51)],[Fmax2(31)/pi/lambda(31)^2,Fmax2(41)/pi/lambda(41)^2,Fmax2(51)/pi/lambda(51)^2],125,c(6:8,:),'Marker','x','LineWidth',2)

axis on
ax = gca;
box on
ax.XLim = [1e-3,1e2];
ax.XScale = "log";
ax.YLim = [0.9,2.1];
set(gca,'FontSize',figure_number_fontsize)
xlabel('$\lambda_E=\sqrt{\Gamma}\,\mathcal{R}_s^2$','interpreter','latex')
ylabel('$\mathcal{F}_c=F_c/(\pi\gamma{R_s})$','interpreter','latex')
xticks(logspace(-3,2,6))
xticklabels({'10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'})
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';

hold off


figure(2)

figure_height = 8*0.75;

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

load('F_H_data_Bradley.mat')

c = parula(10);
Index = zeros(2,length(lambda));
one_line = zeros(1,length(lambda));

plot(H+1,F_Bradley,'color','k','LineWidth',line_width)

for i=1:length(lambda)

    A = isfinite(F_1(:,i));
    B = isfinite(F_2(:,i));

    if numel(find(~A))
        Index(1,i) = length(H)-numel(find(isfinite(F_1(:,i))));
    else
        dF_1 = abs(F_1(2:end,i)-F_1(1:end-1,i))/(pi*lambda_E(i)^2);
        if any(dF_1>0.1)
            AA = find(dF_1>0.1);
            Index(1,i) = AA(1);
        else
            one_line(i) = 1;
        end
    end

    if numel(find(~B))
        Index(2,i) = numel(find(isfinite(F_2(:,i))));
    else
        if one_line(i) == 1
            continue
        else
            dF_2 = abs(F_2(2:end,i)-F_2(1:end-1,i))/(pi*lambda_E(i)^2);
            BB = find(dF_2>0.1);
            Index(2,i) = BB(1);
        end
    end
end

for i=1:length(lambda)
    if one_line(i) == 1
        plot(H+1,F_1(:,i)/(pi*lambda_E(i)^2),'color',c(2+i,:),'LineWidth',line_width)
    else
        plot(H(Index(1,i)+1:end)+1,F_1(Index(1,i)+1:end,i)/(pi*lambda_E(i)^2),'color',c(2+i,:),'linewidth',line_width)
    end
    text(1,0,['$\lambda_E=$',num2str(lambda_E(i))],'Interpreter','latex','color',c(2+i,:))
end

for i=1:length(lambda)
    if one_line(i) == 1
        continue
    else
        plot(H(1:Index(2,i))+1,F_2(1:Index(2,i),i)/pi/lambda_E(i)^2,'color',c(2+i,:),'linewidth',line_width)
    end
end

for i=1:length(lambda)
    if one_line(i) == 1
        continue
    else
        line([H(Index(1,i))+1,H(Index(1,i)+1)+1],[F_2(Index(1,i),i)/pi/lambda_E(i)^2,F_1(Index(1,i)+1,i)/pi/lambda_E(i)^2],'color',c(2+i,:),'linestyle','--')
        line([H(Index(2,i))+1,H(Index(2,i)+1)+1],[F_2(Index(2,i),i)/pi/lambda_E(i)^2,F_1(Index(2,i)+1,i)/pi/lambda_E(i)^2],'color',c(2+i,:),'linestyle','--')
    end
end

axis on
ax = gca;
ax.TickLength = [0.02 0.02];
box on
%ax.XLim = [0.1,100];
ax.XScale = 'log';
ax.YLim = [0,2.5];
set(gca,'FontSize',figure_number_fontsize)
xlabel('$\mathit{\Delta}=\delta/\delta_{\mathrm{micro}}$','interpreter','latex')
ylabel('$\mathcal{F}=F/(\pi\gamma{R_s})$','interpreter','latex')
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';

hold off


figure(3)

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

Position_set = [21.76923076923077,1.164444444444444,1.4e-14;
                9.653846153846157,1.111111111111111,1.4e-14;
                -3.230769230769234 1.084444444444444 1.4e-14];

hold on

load('F_H_data_JKR_lambdaE_1.mat')
Rd_collect = zeros(n,3);
RF_collect = zeros(n,3);
Rd_collect(:,1) = Rd(:,1);
RF_collect(:,1) = RF1(:,1);

load('F_H_data_JKR_lambdaE_10.mat')
Rd_collect(:,2) = Rd(:,1);
RF_collect(:,2) = RF1(:,1);

load('F_H_data_JKR_lambdaE_100.mat')
Rd_collect(:,3) = Rd(:,1);
RF_collect(:,3) = RF1(:,1);

load('F_H_data_approach_JKR.mat')

leg = {' Numerical results',' JKR limit'};
c = parula(10);
Index = zeros(2,length(lambda));
one_line = zeros(1,length(lambda));

for i=1:length(lambda)

    B = isfinite(F_2(:,i));

    if numel(find(~B))
        Index(2,i) = numel(find(isfinite(F_2(:,i))));
    else
        if one_line(i) == 1
            continue
        else 
            dF_2 = abs(F_2(2:end,i)-F_2(1:end-1,i))/(pi*lambda_E(i)^2);
            BB = find(dF_2>0.1);
            Index(2,i) = BB(1);
        end
    end
end

for i=1:length(lambda)

    plot(H(1:Index(2,i))/lambda_E(i),F_2(1:Index(2,i),i)/pi/lambda_E(i)^2,'color',c(5+i,:),'linewidth',line_width)

    plot(Rd_collect(1:25:201,i),RF_collect(1:25:201,i),'color',c(5+i,:),'LineWidth',line_width,'LineStyle','-.')
    plot(Rd_collect(201:15:261,i),RF_collect(201:15:261,i),'color',c(5+i,:),'LineWidth',line_width,'LineStyle','-.')
    plot(Rd_collect(261:5:296,i),RF_collect(261:5:296,i),'color',c(5+i,:),'LineWidth',line_width,'LineStyle','-.')
    plot(Rd_collect(296:300,i),RF_collect(296:300,i),'color',c(5+i,:),'LineWidth',line_width,'LineStyle','-.')

    text(1,0,['$\lambda_E=$',num2str(lambda_E(i))],'Interpreter','latex','color',c(5+i,:),'Position',Position_set(i,:))
end

axis on
ax = gca;
box on
ax.TickLength = [0.02 0.02];
ax.XLim = [-5,40];
ax.YLim = [-0.1,1.5];
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
set(gca,'FontSize',figure_number_fontsize)
xlabel('$\mathit{\bar{\Delta}}=\delta/\delta_{\mathrm{macro}}$','interpreter','latex')
ylabel('$\mathcal{F}=F/(\pi\gamma{R_s})$','interpreter','latex')

%{
    If the file is too large and causes the exported EMF file to appear
    blurry (the file occupy very little memory), try using the following
    statement:
%}

% exportgraphics(gcf, 'figure3_new.emf', 'ContentType', 'vector')

hold off



