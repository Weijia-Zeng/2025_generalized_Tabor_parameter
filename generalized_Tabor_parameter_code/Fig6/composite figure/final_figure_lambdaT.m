clear; clc

figure(1)

figure_width = 8;
figure_height = 12;
figure_number_fontsize = 10;
line_width = 1.5;

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

filename = 'lambda_E_1000_Rf_100_Rs_20';
load(filename)

plot(lambda_T,Fmax/pi./lambda_E.^2,'linewidth',1.5,'Color','k')

% xl = xline(lambda_E,'--k',['$\lambda_E=$',num2str(lambda_E)],'LineWidth',1,'Interpreter','latex','FontSize',12);
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'left';
% xl.LabelOrientation = 'horizontal';

line([1e-3,1e0],[2,2],'linewidth',2.5,'linestyle','--','color','k')
line([1e0,1e3],[1,1],'linewidth',2.5,'linestyle','--','color','k')

c = parula(9);

scatter([lambda(11),lambda(21),lambda(31)],[Fmax(11)/pi/lambda_E^2,Fmax(21)/pi/lambda_E^2,Fmax(31)/pi/lambda_E^2],125,c(3:5,:),'Marker','x','LineWidth',2)
scatter([lambda(41),lambda(51)],[Fmax(41)/pi/lambda_E^2,Fmax(51)/pi/lambda_E^2],125,c(7:8,:),'Marker','x','LineWidth',2)

axis on
ax = gca;
box on
ax.XLim = [1e-3,1e3];
ax.XScale = "log";
ax.YLim = [0.9,2.1];
set(gca,'FontSize',figure_number_fontsize)
xlabel('$\lambda_T=\Gamma \mathcal{R}_s^2/\mathcal{T}_{\mathrm{pre}}$','interpreter','latex')
ylabel('$\mathcal{F}_c=F_c/(\pi\gamma{R_s})$','interpreter','latex')
xticks(logspace(-3,3,4))
xticklabels({'10^{-3}','10^{-1}','10^{1}','10^{3}'})
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';

hold off


figure(2)

figure_height = figure_width*0.75;
set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

load('F_H_data_Bradley.mat')

c = parula(9);
Index = zeros(2,length(lambda));
one_line = zeros(1,length(lambda));

plot(H+1,F_Bradley,'color','k','LineWidth',line_width)

for i=1:length(lambda) 

    A = isfinite(F_1(:,i));
    B = isfinite(F_2(:,i));

    if numel(find(~A))
        Index(1,i) = length(H)-numel(find(isfinite(F_1(:,i))));
    else
        dF_1 = abs(F_1(2:end,i)-F_1(1:end-1,i))/(pi*lambda_E^2);
        if any(dF_1>0.18)
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
            dF_2 = abs(F_2(2:end,i)-F_2(1:end-1,i))/(pi*lambda_E^2);
            BB = find(dF_2>0.1);
            Index(2,i) = BB(1);
        end
    end
end

for i=1:length(lambda)
    if one_line(i) == 1
        plot(H+1,F_1(:,i)/(pi*lambda_E^2),'color',c(2+i,:),'LineWidth',line_width)
    else
        plot(H(Index(1,i)+1:end)+1,F_1(Index(1,i)+1:end,i)/(pi*lambda_E^2),'color',c(2+i,:),'linewidth',line_width)
    end
    text(1,0,['$\lambda_T=$',num2str(lambda_T(i))],'interpreter','latex','Color',c(2+i,:))
end

for i=1:length(lambda)
    if one_line(i) == 1
        continue
    else
        plot(H(1:Index(2,i))+1,F_2(1:Index(2,i),i)/pi/lambda_E^2,'color',c(2+i,:),'linewidth',line_width)
    end
end

for i=1:length(lambda)
    if one_line(i) == 1
        continue
    else
        line([H(Index(1,i))+1,H(Index(1,i)+1)+1],[F_2(Index(1,i),i)/pi/lambda_E^2,F_1(Index(1,i)+1,i)/pi/lambda_E^2],'color',c(2+i,:),'linestyle','--')
        line([H(Index(2,i))+1,H(Index(2,i)+1)+1],[F_2(Index(2,i),i)/pi/lambda_E^2,F_1(Index(2,i)+1,i)/pi/lambda_E^2],'color',c(2+i,:),'linestyle','--')
    end
end

axis on
ax = gca;
ax.TickLength = [0.02 0.02];
box on
%ax.XLim = [0.1,100];
ax.XScale = 'log';
ax.YLim = [0,2.5];
set(gca,'Fontsize',figure_number_fontsize)
xlabel('$\mathit{\Delta}=\delta/\delta_{\mathrm{micro}}$','interpreter','latex')
ylabel('$\mathcal{F}=F/(\pi\gamma{R_s})$','interpreter','latex')
% xticks([1 10])
% xticklabels({'10^{0}','10^{1}'})
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';


figure(3)

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

load('F_H_data_JKR_lambdaT_1.mat')
Rd_collect = zeros(n,3);
RF_collect = zeros(n,3);
Rd_collect(:,1) = Rd(:,1);
RF_collect(:,1) = RF1(:,1);

load('F_H_data_JKR_lambdaT_10.mat')
Rd_collect(:,2) = Rd(:,1);
RF_collect(:,2) = RF1(:,1);

load('F_H_data_JKR_lambdaT_100.mat')
Rd_collect(:,3) = Rd(:,1);
RF_collect(:,3) = RF1(:,1);

load('F_H_data_approach_JKR.mat')

leg = {' Numerical results',' JKR limit'};
c = parula(9);
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
            dF_2 = abs(F_2(2:end,i)-F_2(1:end-1,i))/(pi*lambda_E^2);
            BB = find(dF_2>0.1);
            Index(2,i) = BB(1);
        end
    end
end

% plot(H(1:Index(2,1))/lambda_E(1),F_2(1:Index(2,1),1)/pi/lambda_E(1)^2,'color','k','linewidth',2.5);
% plot(Rd_collect(:,1)/lambda_E(1),RF_collect(:,1),'color','k','LineWidth',2,'LineStyle','--');

for i=1:length(lambda)

    if i == 1
        plot(H(1:Index(2,i))/lambda(i),F_2(1:Index(2,i),i)/pi/lambda_E^2,'color',c(4+i,:),'linewidth',line_width)
        plot(Rd_collect(:,i),RF_collect(:,i),'color',c(4+i,:),'LineWidth',line_width,'LineStyle','--')
        text(1,0.8,['$\lambda_T=$',num2str(lambda_T(i))],'interpreter','latex','Color',c(4+i,:))
    else
        plot(H(1:Index(2,i))/lambda(i),F_2(1:Index(2,i),i)/pi/lambda_E^2,'color',c(5+i,:),'linewidth',line_width)
        plot(Rd_collect(:,i),RF_collect(:,i),'color',c(5+i,:),'LineWidth',line_width,'LineStyle','--')
        text(1,0.8,['$\lambda_T=$',num2str(lambda_T(i))],'interpreter','latex','Color',c(5+i,:))
    end
    
end

axis on
ax = gca; 
ax.TickLength = [0.02 0.02];
box on
ax.XLim = [0,3];
% ax.XScale = "log";
ax.YLim = [0.5,1.5];
set(gca,'Fontsize',figure_number_fontsize)
xlabel('$\mathit{\bar{\Delta}}=\delta/\delta_{\mathrm{macro}}$','interpreter','latex')
ylabel('$\mathcal{F}=F/(\pi\gamma{R_s})$','interpreter','latex')
ax.XMinorTick = 'off';
ax.YMinorTick = 'off'; 



