clear; clc

load('F_H_data_approach_JKR.mat')

Rs_exact = 400;   % Unit: nm
z0_exact = 1;     % Unit: nm

R_link_exact = [0 Rmesh] * sqrt(Rs_exact*z0_exact);
R_link = [0 Rmesh];
gap_link = zeros(length(lambda_E),n);

w_position = zeros(1,length(lambda_E));
van_der_Waal_position = zeros(1,length(lambda_E));
w_max = zeros(1,length(lambda_E));
van_der_Waal_max = zeros(1,length(lambda_E));

default_color = {'#0072BD','#D95319','#EDB120'};
c = parula(10);
leg = cell(1,length(lambda_E));

%% membrane deformation
h1 = figure(1);

figure_width = 7.5;
figure_height = figure_width*0.75;
figure_number_fontsize = 10;
line_width = 1.5;

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

for i = 1:length(lambda_E)
    w_link_exact = [(4*w_pulloff(1,i)-w_pulloff(2,i))/3 w_pulloff(:,i)']*z0_exact;

    [w_max_exact(i),w_position(i)] = max(w_link_exact);

    R_link_exact_2 = [-flip(Rmesh) 0 Rmesh] * sqrt(Rs_exact*z0_exact);
    w_link_exact_2 = [flip( w_pulloff(:,i)' ) (4*w_pulloff(1,i)-w_pulloff(2,i))/3 w_pulloff(:,i)']*z0_exact;

    plot(R_link_exact,w_link_exact,'LineWidth',line_width,'Color',c(5+i,:))
    rectangle('Position',[-Rs_exact,( H_pulloff(i)+1 )*z0_exact,2*Rs_exact,2*Rs_exact],'Curvature',[1,1],'FaceColor','none','EdgeColor',c(5+i,:),'LineWidth',line_width)
    
    leg{i} = ['$\lambda_E=$',num2str(lambda_E(i))];
end

ax = gca;
axis equal
set(gca,'fontsize',figure_number_fontsize)
xlabel('$r/{\rm nm}$','Interpreter','latex')
ylabel('$w/{\rm nm}$','Interpreter','latex')
% legend(leg,'Interpreter','latex','Location','best','Box','off')
box on
ax.XLim = [0,2000];

zp = BaseZoom([950,500,750,750],[160,500,120,120]);
zp.run
zp.subAxes.Color = 'white';
zp.subAxes.YAxisLocation = "right";

set(gca,'looseInset',[0 0 0 0.03])

% close(h1)

%% surface gap
h2 = figure(2);

figure_width = 7.5;
figure_height = figure_width*0.75;
set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

for i = 1:length(lambda_E)
    gap = H_pulloff(i)+1+Rmesh.^2/2-w_pulloff(:,i)';
    gap_link(i,:) = [H_pulloff(i)+1-(4*w_pulloff(1,i)-w_pulloff(2,i))/3 gap];

    plot(R_link_exact,gap_link(i,:),'color',c(5+i,:),'LineWidth',line_width,'LineStyle','-')
end

for i = 1:length(lambda_E)
    scatter(R_link_exact(w_position(i)),gap_link(i,w_position(i)),[],c(5+i,:))
end

yl = yline(3^(1/6),'LineWidth',line_width,'LineStyle','--');
yl.Label = '''contact line''';
yl.LabelVerticalAlignment = 'top';
yl.LabelHorizontalAlignment = 'center';
yl.LabelOrientation = 'horizontal';

ax = gca;

ax.XLim = [0,250];
ax.YLim = [0.9,1.5];

box on
set(gca,'FontSize',figure_number_fontsize)
xlabel('$r/{\rm nm}$','Interpreter','latex')
ylabel('$g/{\rm nm}$','Interpreter','latex')

% close(h2)

%% van der Waals force
h3 = figure(3);

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

for i = 1:length(lambda_E)

    van_der_Waal = 8/3*lambda_E(i)^2/Rs^2*(gap_link(i,:).^(-3)-gap_link(i,:).^(-9));
    [van_der_Waal_max(i),van_der_Waal_position(i)] = max(van_der_Waal);
    plot(R_link_exact,van_der_Waal/van_der_Waal_max(i),'Color',c(5+i,:),'LineWidth',line_width)

end

% for i = 1:length(lambda_E)
% 
%     xline(R_link(van_der_Waal_position(i)),'LineStyle','--');
% 
% end

ax = gca;
ax.XLim = [0,250];
ax.YLim = [-0.8,1];

box on
set(gca,'FontSize',figure_number_fontsize)
xlabel('$r/{\rm nm}$','Interpreter','latex')
ylabel('$q(r)/q_{\rm max}$','Interpreter','latex')

close(h3)

%% adhesion density
h4 = figure(4);

set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

hold on

for i = 1:length(lambda_E)

        U = 8/3*lambda_E(i)^2/Rs^2*(-1/2*gap_link(i,:).^(-2)+1/8*gap_link(i,:).^(-8));
        plot(R_link_exact,U/(-lambda_E(i)^2/Rs^2),'Color',c(5+i,:),'LineWidth',line_width)

end

D = 3^(1/6);
U_mark = -8/3*( -1/(2*D^2) + 1/(8*D^8) );

yl = yline(U_mark,'LineWidth',line_width,'LineStyle','--');
yl.Label = '''contact line''';
yl.LabelVerticalAlignment = 'bottom';
yl.LabelHorizontalAlignment = 'center';
yl.LabelOrientation = 'horizontal';

ax = gca;
ax.XLim = [0,250];
ax.YLim = [0,1.05];

box on
set(gca,'FontSize',figure_number_fontsize)
xlabel('$r/{\rm nm}$','Interpreter','latex')
ylabel('$U_{\mathrm{LJ}}/\gamma$','Interpreter','latex')

% close(h4)

%% ratio of energy
h5 = figure(5);

figure_height = figure_width*0.75;
set(gcf,"Units","centimeters")
set(gcf,'Position',[10 10 figure_width figure_height])

U_adhesion = zeros(1,length(lambda_E));
U_JKR = zeros(1,length(lambda_E));

dR = Rf/n;

hold on

for i = 1:length(lambda_E)

    % for j = 1:van_der_Waal_position(i)
    for j = 1:n

        D = gap_link(i,j);
        U = 8/3*lambda_E(i)^2/Rs^2*( -1/(2*D^2) + 1/(8*D^8) );

        % if j == 1 || j == van_der_Waal_position(i)
        if j == 1 || j == n
            U_adhesion(i) = U_adhesion(i) + 1/2*U*R_link(j)*2*pi*dR;
        else
            U_adhesion(i) = U_adhesion(i) + U*R_link(j)*2*pi*dR;
        end

    end

    U_JKR(i) = -lambda_E(i)^2/Rs^2*pi*R_link(van_der_Waal_position(i))^2;

end

for i = 1:length(lambda_E)

    Bar = bar(i,U_adhesion(i)./U_JKR(i));
    Bar.FaceColor = c(5+i,:);
    Bar.EdgeColor = 'none';
    text(Bar.XEndPoints,Bar.YEndPoints,num2str(U_adhesion(i)./U_JKR(i)),'HorizontalAlignment','center','VerticalAlignment','bottom')

end

ax = gca;
ax.XTick = [1,2,3];
ax.XTickLabel = {'1','10','100'};

box on
set(gca,'FontSize',figure_number_fontsize)
xlabel('$\lambda_E=\sqrt{\Gamma}\,\mathcal{R}_s^2$','Interpreter','latex')
ylabel('$U_{\rm ad}/(\pi\gamma a^2)$','Interpreter','latex')
% set(gca,'looseInset',[0 0 0.03 0.05])

% close(h5)



