clear; clc

%%%%%%%%%%%%%%%%  dimensionless parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%
Rf_total = [20 100];
Rs_total = [10 20];
T_total = [0 0.01 0.1];
P_total = [0 1e-6];

%% numerical solutions figure
figure(1)

figure_width = 8;
figure_height = 8*0.75;
figure_number_fontsize = 10;
line_width = 1.5;

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

colormap parula

for ii = 1:length(Rf_total)

    for jj = 1:length(Rs_total)

        for kk = 1:length(T_total)

            for ll = 1:length(P_total)

                load(['Rf_',num2str(Rf_total(ii)),'_Rs_',num2str(Rs_total(jj)),'_T_',num2str(T_total(kk)),'_P_',num2str(P_total(ll)),'_new.mat'])

                if ii == 1
                    sz = 10;    % Rf = 20
                else
                    sz = 20;    % Rf = 100
                end

                if jj == 1
                    s = scatter(lambda,Fmax,sz,log10(Gamma)); % Rs = 10
                else
                    s = scatter(lambda,Fmax,sz,log10(Gamma),'filled'); % Rs = 20
                end

                if ll == 2
                    s.MarkerFaceAlpha = 0.5; % pressure exists
                end

                if kk == 1
                    s.Marker = 'o'; % T = 0
                elseif kk == 2
                    s.Marker = '^'; % T = 1
                else
                    s.Marker = 'square';  % T = 10
                end

            end
        end
    end
end

axis on
ax = gca;
box on
% ax.XGrid = 'on';
% ax.YGrid = 'on';
% ax.MinorGridAlpha = 0.1;
% ax.GridAlpha = 0.1;
% ax.TickLength = [0.02 0.02];
ax.XLim = [1e-3,10^(1.5)];
% ax.XLim = [1e-3,1e2];
ax.XTick = logspace(-3,2,6);
ax.XScale = "log";
ax.YLim = [0.85,2.15];
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
set(gca,'FontSize',figure_number_fontsize)
xlabel('$\lambda$','interpreter','latex')
ylabel('$\mathcal{F}_c=F_c/(\pi\gamma{R_s})$','interpreter','latex')

% used to modify the width of the colorbar
c = colorbar('Location','eastoutside');
c.Label.Interpreter = 'latex';
c.Label.String = '$\log_{10} \Gamma$';
c.FontSize = figure_number_fontsize;
ax.Position(3) = 0.645;
axposition = ax.Position;
cposition = c.Position;
cposition(3) = 0.5*cposition(3);
c.Position = cposition;
ax.Position = axposition;

hold off

%% compared with experimental data
figure(2)

figure_width = 8;
figure_height = 8*0.75;
figure_number_fontsize = 10;
line_width = 1.5;

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

load('Exp_data.mat')    % load the previous experimental data

ga = [0.049 0.049 0.060 0.060 0.060 0.060 0.050 0 0.055 0.041 0.041];  % gamma from previous work

c = parula(6);
c_index = [1 2 3 3 3 3 4 0 5 5 6];

for i = [1 2 3 4 5 6 7 9 10 11]     % remove the eighth set of data
    lambda = ( -Exp_data(i).tpre+sqrt(Exp_data(i).tpre.^2+4*1000*Exp_data(i).t*ga(i)) )/(2*1000*Exp_data(i).t)*Exp_data(i).Radius/1;
    Fout = -Exp_data(i).Fout/pi/ga(i)./Exp_data(i).Radius;

    switch i
        case 1
            h1 = scatter(lambda,Fout,[],c(c_index(i),:),'filled');
        case 2
            h2 = scatter(lambda,Fout,[],c(c_index(i),:),'filled');
        case 3
            h3 = scatter(lambda,Fout,[],c(c_index(i),:),'filled');
        case 7
            h4 = scatter(lambda,Fout,[],c(c_index(i),:),'filled');
        case 9
            h5 = scatter(lambda,Fout,[],c(c_index(i),:),'filled');
        case 11
            h6 = scatter(lambda,Fout,[],c(c_index(i),:),'filled');
        otherwise
            scatter(lambda,Fout,[],c(c_index(i),:),'filled');
    end
end

colormap parula

% lower bound
load('Rf_100_Rs_20_T_0_P_0_new.mat')
lambda1 = lambda;
Fmax_low = Fmax;
plot(lambda,Fmax,'color',[168 208 219]/258)

% upper bound
load('Rf_20_Rs_20_T_0.1_P_1e-06_new.mat')
lambda2 = fliplr(lambda);
Fmax_high = fliplr(Fmax);
plot(lambda,Fmax,'color',[168 208 219]/258)

lambda_link = [lambda1 lambda2];
Fmax_link = [Fmax_low Fmax_high];
patch(lambda_link,Fmax_link,[168 208 219]/258,'FaceAlpha',0.3,'EdgeColor','none')


axis on
ax = gca;
box on
% ax.XGrid = 'on';
% ax.YGrid = 'on';
% ax.MinorGridAlpha = 0.1;
% ax.GridAlpha = 0.1;
% ax.TickLength = [0.02 0.02];
ax.XLim = [1e-3,10^(1.5)];
% ax.XLim = [1e-3,1e2];
ax.XTick = logspace(-3,2,6);
ax.XScale = "log";
ax.YLim = [0.85,2.15];
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
set(gca,'FontSize',figure_number_fontsize)
xlabel('$\lambda$','interpreter','latex')
ylabel('$\mathcal{F}_c=F_c/(\pi\gamma{R_s})$','interpreter','latex')

%% legned
Leg = legend([h1 h2 h3 h4 h5 h6],{'$2{\rm L},R_f \approx 1800{\rm nm},Rs \approx 40{\rm nm}$','$2{\rm L},R_f \approx 2800{\rm nm},Rs \approx 40{\rm nm}$', ...
    '$4{\rm L},R_f \approx 1800/2800{\rm nm},Rs \approx 90{\rm nm}$','$2{\rm L},R_f \approx 2800{\rm nm},Rs \approx 90{\rm nm}$' ...
    '$4{\rm L},R_f \approx 1800/2800{\rm nm},Rs \approx 900{\rm nm}$','$4{\rm L},R_f \approx 1800{\rm nm},Rs \approx 900{\rm nm}$'},'interpreter','latex','Location','best','box','off' ...
    ,'FontSize',7);

pause(1e-6)     % This sentence is indispensable
for i = 1:length(Leg.EntryContainer.NodeChildren)
    Leg.EntryContainer.NodeChildren(i).Children(1).Transform.Children.Children.VertexData = single([0.9;0.5;0]);  % horizontal position, vertical position and in the third direction for the marker
end

print('experiment_data','-dmeta')

hold off
