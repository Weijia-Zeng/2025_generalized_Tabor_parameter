lambdaE = 1;
Rf_filename = [20 100 500];
color_set = {'#E49273','#A8D0DB','#A37A74'};
line_set = cell(1,length(Rf_filename));
Leg = cell(1,length(Rf_filename));

figure(1)

hold on

figure_width = 8;
figure_height = 8*0.75;
figure_number_fontsize = 8;
line_width = 1.5;

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

fill([1e-3,lambdaE,lambdaE,1e-3],[0.85,0.85,2.15,2.15],[255 242 204]/255,'Edgecolor','none','FaceAlpha',0.6)
fill([lambdaE,1e3,1e3,lambdaE],[0.85,0.85,2.15,2.15],[222 235 247]/255,'Edgecolor','none','FaceAlpha',0.6)

fill([1e-3,lambdaE,lambdaE,1e-3],[1,1,2,2],[255 242 204]/255,'Edgecolor','none','FaceAlpha',0.6)
fill([lambdaE,1e3,1e3,lambdaE],[1,1,2,2],[222 235 247]/255,'Edgecolor','none','FaceAlpha',0.6)

for ii = 1:length(Rf_filename)
    filename = ['lambda_E_',num2str(lambdaE),'_Rf_',num2str(Rf_filename(ii)),'_Rs_20.mat'];
    load(filename)
    line_set{ii} = plot(lambda_p,Fmax/pi/lambda_E^2,'Color',color_set{ii},'LineWidth',line_width);
    Leg{ii} = ['${\mathcal R}_f=$',num2str(Rf_filename(ii))];
    % plot(lambda_p,Fmax/pi/lambda_E^2.*(1+2*w_center/Rf^2),'Color','#A8D0DB')
end

F_Bradley_max = zeros(1,length(lambda));

for ii = 1:length(Rf_filename)

    filename = ['lambda_E_',num2str(lambdaE),'_Rf_',num2str(Rf_filename(ii)),'_Rs_20.mat'];
    load(filename)

    for i = 1:length(lambda)

        F_Bradley = zeros(1,length(H));
    
        for j = 1:length(H)
    
            for k=1:n-1
                F_Bradley(j) = F_Bradley(j)+16/3*pi*lambda_E^2*Rmesh(k)*dR*((w_center(i)+H(j)+Rmesh(k)^2/2-w_p(k,i)+1)^(-3)-(w_center(i)+H(j)+Rmesh(k)^2/2-w_p(k,i)+1)^(-9));
            end
            F_Bradley(j) = F_Bradley(j)+16/3*pi*lambda_E^2*Rf*dR*((w_center(i)+H(j)+Rf^2/2-wn+1)^(-3)-(w_center(i)+H(j)+Rf^2/2-wn+1)^(-9));
    
            if j>1 && F_Bradley(j)<F_Bradley(j-1)
    
                F_Bradley_max(i) = F_Bradley(j-1);
                break
    
            end
        end
    end

    plot(lambda_p,F_Bradley_max/pi/lambda_E^2,'Color',color_set{ii},'LineWidth',line_width,'Linestyle','-.')
end

% xl = xline(lambda_E,'--','FontSize',figure_number_fontsize);
% xl.Interpreter = "latex";
% xl.Label = ['$\lambda_E=$',num2str(lambda_E)];
% xl.LabelVerticalAlignment = 'bottom';
% xl.LabelHorizontalAlignment = 'center';
% xl.LabelOrientation = 'horizontal';

filename = ['lambda_E_',num2str(lambdaE),'_Rf_',num2str(Rf_filename(2)),'_Rs_20.mat'];
load(filename)

% yl = yline(Fmax(end)/pi/lambda_E^2,'--k',num2str(Fmax(end)/pi/lambda_E^2),'Interpreter','latex','FontSize',figure_number_fontsize,'LineWidth',line_width);
% yl.LabelVerticalAlignment = 'top';
% yl.LabelHorizontalAlignment = 'left';
% yl.LabelOrientation = 'horizontal';

line([1e-3,10^(-1.5)],[2,2],'linestyle','--','color','k','Linewidth',line_width)
line([10^(1.5),1e3],[1,1],'linestyle','--','color','k','Linewidth',line_width)

% text(1e-2,1,['$\lambda_E=$',num2str(lambdaE)],'interpreter','latex','Position',[0.325702087916432,0.926553673061442,0])
% text(1e-2,1,'$\lambda_p ~ \rm dom$','interpreter','latex','Position',[0.00701703876823,2.072316384925849,0])
% text(1e-2,1,'$\lambda_E ~ \rm dom$','interpreter','latex','Position',[11.253356598239787,2.072316384925849,0])

color = parula(10);
load('Rf_100_Rs_20_reference.mat')
scatter(1e3,Fmax(31)/pi/lambda(31)^2,100,color(6,:),'Marker','x','LineWidth',2)

% annotation('textbox',[0.615894039735099,0.277533039647578,0.21523178807947,0.096916299559472],'String',['$\mathcal{R}_f=$',num2str(Rf)],'Interpreter','latex')

% c = [.01 .72 .77
%     .99 .49 .00
%     .68 .42 .89];

% scatter([1e-2,10^(-0.5),10],[Fmax(11)/pi/lambda(11)^2,Fmax(26)/pi/lambda(26)^2,Fmax(41)/pi/lambda(41)^2],75,c,'Marker','x','LineWidth',2)

axis on
box on
ax = gca;
ax.XScale = "log";
ax.XTick = logspace(-3,3,7);
ax.XLim = [1e-3,1e3];
ax.YLim = [0.85 2.15];
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.XMinorGrid = 'off';
ax.YMinorGrid= 'off';
ax.MinorGridAlpha = 0.1;
ax.GridAlpha = 0.1;
ax.TickLength = [0.02 0.02];
set(gca,'FontSize',figure_number_fontsize)
xlabel('$\lambda_p=\Gamma\mathcal{R}_s^{8/3}/(P^{2/3}\mathcal{R}_f^{2/3})$','interpreter','latex')
ylabel('$\mathcal{F}_c=F_c/(\pi\gamma{R_s})$','interpreter','latex')
% title(['$\mathcal{R}_f=$',num2str(Rf),'$\quad \mathcal{R}_s=$',num2str(Rs),'$\quad \lambda_E=$',num2str(lambda_E),'$\quad \lambda_T=\infty$'],'Interpreter','latex')
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',line_width)
legend([line_set{1:end}],Leg,'Interpreter','latex',"Box","off",'Position',[0.57276079073087,0.536588508374271,0.291929531851083,0.174273122963926])

hold off