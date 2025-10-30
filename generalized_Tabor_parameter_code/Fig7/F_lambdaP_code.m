clear; clc

n = 5000;
start = -0.3;
cutoff = 1000;
H = linspace(start,cutoff,(cutoff-start)*100+1);
v = 0.165;

Rf = 20;
Rs = 20;
lambda_E = 0.1;

start = -3;
cutoff = 3;
lambda_p = logspace(start,cutoff,ceil( (cutoff-start)*10 )+1);
p = lambda_E^3./lambda_p.^(1.5)/Rf;   % dimensionless pressure
% lambda_T = NaN;
% T = lambda_E^2/lambda_T;
% lambda = -lambda_E^2/2*(1./lambda_p+1/lambda_T)+sqrt( (lambda_E^2/2*(1./lambda_p+1/lambda_T)).^2+lambda_E^2 );
T = 0;
lambda = -lambda_E^2/2*(1./lambda_p)+sqrt( (lambda_E^2/2*(1./lambda_p)).^2+lambda_E^2 );

W = 0.1*Rf*(1-((1:n-1)'./n).^2);        
PHI = 0.001*ones(n-1,1);
beta = 1;

w_p = zeros(n-1,length(lambda_p));
phi_p = zeros(n-1,length(lambda_p));
w_center = zeros(1,length(lambda_p));


%%%%%%%%%%%%%%%%  start to solve  %%%%%%%%%%%%%%%%%%%

tic
f = zeros(n-1,1);
g = zeros(n-1,1);
e = zeros(2*n-2,2*n-2);
b = zeros(2*n-2,1);
m = zeros(2*n-2,1);

w_pulloff = zeros(n-1,length(lambda));
phi_pulloff = zeros(n-1,length(lambda));
H_pulloff = zeros(1,length(lambda));

dR = Rf/n;
Rmesh = (1:n-1).*dR;
     
F = zeros(length(H),length(lambda));
Fmax = zeros(1,length(lambda));


%%%%%%%%%%%%%%%%% pressurize bubble shape %%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(p)

    R = ones(2*n-2,1);
    w = W;
    phi = PHI;

    while(norm(R)>1e-9)

        % boundary condition
        w0 = w(1);
        wn = 0;
        phi0 = phi(1)-v*dR*phi(1)/Rmesh(1);
        phin = ( phi(n-1)+(1-v)*dR*T ) / ( 1-v*dR/Rf );

        for k=2:(n-2) 
            f(k) = phi(k)/Rmesh(k) * (w(k+1)-2*w(k)+w(k-1))/dR^2 + 1/Rmesh(k) * (phi(k+1)-phi(k-1))/(2*dR) * (w(k+1)-w(k-1))/(2*dR) + p(i);
            g(k) = (phi(k+1)-phi(k-1))/(2*dR) - phi(k)/Rmesh(k) + Rmesh(k) * (phi(k+1)-2*phi(k)+phi(k-1))/dR^2 + 1/2 * ( (w(k+1)-w(k-1))/(2*dR) )^2;

            e(k,k-1) = phi(k)/Rmesh(k) * 1/dR^2 - 1/Rmesh(k) * (phi(k+1)-phi(k-1))/(2*dR) * 1/(2*dR);
            e(k,k) = phi(k)/Rmesh(k) * (-2)/dR^2;
            e(k,k+1) = phi(k)/Rmesh(k) * 1/dR^2 + 1/Rmesh(k) * (phi(k+1)-phi(k-1))/(2*dR) * 1/(2*dR);
            e(k,k+n-2) = -1/Rmesh(k) * 1/(2*dR) * (w(k+1)-w(k-1))/(2*dR);
            e(k,k+n-1) = 1/Rmesh(k) * (w(k+1)-2*w(k)+w(k-1))/dR^2;
            e(k,k+n) = 1/Rmesh(k) * 1/(2*dR) * (w(k+1)-w(k-1))/(2*dR);
      
            e(k+n-1,k-1) = (w(k+1)-w(k-1))/(2*dR) * (-1)/(2*dR);
            e(k+n-1,k) = 0;
            e(k+n-1,k+1) = (w(k+1)-w(k-1))/(2*dR) * 1/(2*dR);
            e(k+n-1,k+n-2) = -1/(2*dR) + Rmesh(k) * 1/dR^2;
            e(k+n-1,k+n-1) = -1/Rmesh(k) + Rmesh(k) * (-2)/dR^2;
            e(k+n-1,k+n) = 1/(2*dR) + Rmesh(k) * 1/dR^2;
        end

        f(1) = phi(1)/Rmesh(1) * (w(2)-2*w(1)+w0)/dR^2 + p(i);
        g(1) = (phi(1)-phi0)/dR - phi(1)/Rmesh(1) + Rmesh(1) * (phi(2)-2*phi(1)+phi0)/dR^2;

        e(1,1) =  phi(1)/Rmesh(1) * (-2+1)/dR^2;
        e(1,2) = phi(1)/Rmesh(1) * 1/dR^2;
        e(1,n) = 1/Rmesh(1) * (w(2)-2*w(1)+w0)/dR^2;
        e(1,n+1) = 0;

        e(n,1) = 0;
        e(n,2) = 0;
        e(n,n) = (1-1+v*dR/Rmesh(1))/dR - 1/Rmesh(1) + Rmesh(1) * (-2+1-v*dR/Rmesh(1))/dR^2;
        e(n,n+1) = Rmesh(1) * 1/dR^2;

        f(n-1) = phi(n-1)/Rmesh(n-1) * (wn-2*w(n-1)+w(n-2))/dR^2 + 1/Rmesh(n-1) * (phin-phi(n-2))/(2*dR) * (wn-w(n-2))/(2*dR) + p(i);
        g(n-1) = (phin-phi(n-1))/dR - phi(n-1)/Rmesh(n-1) + Rmesh(n-1) * (phin-2*phi(n-1)+phi(n-2))/dR^2 + 1/2 * ( (wn-w(n-1))/dR )^2;

        e(n-1,n-2) = phi(n-1)/Rmesh(n-1) * 1/dR^2 + 1/Rmesh(n-1) * (phin-phi(n-2))/(2*dR) * (-1)/(2*dR);
        e(n-1,n-1) = phi(n-1)/Rmesh(n-1) * (-2)/dR^2;
        e(n-1,2*n-3) = 1/Rmesh(n-1) * (-1)/(2*dR) * (wn-w(n-2))/(2*dR);
        e(n-1,2*n-2) = 1/Rmesh(n-1) * (wn-2*w(n-1)+w(n-2))/dR^2 + 1/Rmesh(n-1) * (1/(1-v*dR/Rf))/(2*dR) * (wn-w(n-2))/(2*dR);

        e(2*n-2,n-2) = 0;
        e(2*n-2,n-1) = w(n-1)/dR^2;
        e(2*n-2,2*n-3) = Rmesh(n-1) * 1/dR^2;
        e(2*n-2,2*n-2) = (1/(1-v*dR/Rf)-1)/dR - 1/Rmesh(n-1) + Rmesh(n-1) * (1/(1-v*dR/Rf)-2)/dR^2;
   
        m(1:n-1) = w(1:n-1);
        m(n:2*n-2) = phi(1:n-1);
   
        b(1:n-1) = f(1:n-1);
        b(n:2*n-2) = g(1:n-1);
   
        e = sparse(e);
        delta = -e\b;
   
        R = delta./m;
   
        m = m + beta * delta;
   
        w(1:n-1) = m(1:n-1);
        phi(1:n-1) = m(n:2*n-2);
    end

    w_p(:,i) = w;
    phi_p(:,i) = phi;
    w_center(i) = w_p(1,i);

end

ratio_1 = w_center/Rs/Rf;
ratio_2 = w_center/Rf^2;


%%%%%%%%%%%%%%%%%%%% Lifting the sphere %%%%%%%%%%%%%%%%%%%%%%%
beta = 1;

flag = 0;
iteration_max = 1000;
iteration_middle = iteration_max/2;
warning_index = zeros(1,length(lambda_E));

for i = 1:length(p)

    if i == 1
        w = w_p(:,i);
        phi = phi_p(:,i);
    end
  
    for j = 1:length(H)
      
        R = ones(2*n-2,1);
        iteration = 0;
        if j < flag
            continue
        end
    
        while(norm(R)>1e-9 && iteration ~= iteration_max)
    
            % avoid the membrane penetrateing the sphere surface
            for k=1:n-1
                if w(k) >= w_center(i)+H(j)+0.98+1/2*Rmesh(k)^2
                    w(k) = w_center(i)+H(j)+0.98+1/2*Rmesh(k)^2;
                end
            end

            % boundary conditions
            w0 = w(1);
            wn = 0;
            phi0 = phi(1)-v*dR*phi(1)/Rmesh(1);
            phin = ( phi(n-1)+(1-v)*dR*T ) / ( 1-v*dR/Rf );

            for k=2:(n-2) 
                f(k) = phi(k)/Rmesh(k) * (w(k+1)-2*w(k)+w(k-1))/dR^2 + 1/Rmesh(k) * (phi(k+1)-phi(k-1))/(2*dR) * (w(k+1)-w(k-1))/(2*dR) + p(i) + 8*lambda_E^2/3*((w_center(i)+H(j)+1+Rmesh(k)^2/2-w(k))^(-3)-(w_center(i)+H(j)+1+Rmesh(k)^2/2-w(k))^(-9));
                g(k) = (phi(k+1)-phi(k-1))/(2*dR) - phi(k)/Rmesh(k) + Rmesh(k) * (phi(k+1)-2*phi(k)+phi(k-1))/dR^2 + 1/2 * ( (w(k+1)-w(k-1))/(2*dR) )^2;

                e(k,k-1) = phi(k)/Rmesh(k) * 1/dR^2 - 1/Rmesh(k) * (phi(k+1)-phi(k-1))/(2*dR) * 1/(2*dR);
                e(k,k) = phi(k)/Rmesh(k) * (-2)/dR^2 + 8*lambda_E^2/3*(3*(w_center(i)+H(j)+1+Rmesh(k)^2/2-w(k))^(-4)-9*(w_center(i)+H(j)+1+Rmesh(k)^2/2-w(k))^(-10));
                e(k,k+1) = phi(k)/Rmesh(k) * 1/dR^2 + 1/Rmesh(k) * (phi(k+1)-phi(k-1))/(2*dR) * 1/(2*dR);
                e(k,k+n-2) = -1/Rmesh(k) * 1/(2*dR) * (w(k+1)-w(k-1))/(2*dR);
                e(k,k+n-1) = 1/Rmesh(k) * (w(k+1)-2*w(k)+w(k-1))/dR^2;
                e(k,k+n) = 1/Rmesh(k) * 1/(2*dR) * (w(k+1)-w(k-1))/(2*dR);
      
                e(k+n-1,k-1) = (w(k+1)-w(k-1))/(2*dR) * (-1)/(2*dR);
                e(k+n-1,k) = 0;
                e(k+n-1,k+1) = (w(k+1)-w(k-1))/(2*dR) * 1/(2*dR);
                e(k+n-1,k+n-2) = -1/(2*dR) + Rmesh(k) * 1/dR^2;
                e(k+n-1,k+n-1) = -1/Rmesh(k) + Rmesh(k) * (-2)/dR^2;
                e(k+n-1,k+n) = 1/(2*dR) + Rmesh(k) * 1/dR^2;
            end

            f(1) = phi(1)/Rmesh(1) * (w(2)-2*w(1)+w0)/dR^2 + p(i) + 8*lambda_E^2/3*((w_center(i)+H(j)+1+Rmesh(1)^2/2-w(1))^(-3)-(w_center(i)+H(j)+1+Rmesh(1)^2/2-w(1))^(-9));
            g(1) = (phi(1)-phi0)/dR - phi(1)/Rmesh(1) + Rmesh(1) * (phi(2)-2*phi(1)+phi0)/dR^2;

            e(1,1) =  phi(1)/Rmesh(1) * (-2+1)/dR^2 + 8*lambda_E^2/3*(3*(w_center(i)+H(j)+1+Rmesh(1)^2/2-w(1))^(-4)-9*(w_center(i)+H(j)+1+Rmesh(1)^2/2-w(1))^(-10));
            e(1,2) = phi(1)/Rmesh(1) * 1/dR^2;
            e(1,n) = 1/Rmesh(1) * (w(2)-2*w(1)+w0)/dR^2;
            e(1,n+1) = 0;

            e(n,1) = 0;
            e(n,2) = 0;
            e(n,n) = (1-1+v*dR/Rmesh(1))/dR - 1/Rmesh(1) + Rmesh(1) * (-2+1-v*dR/Rmesh(1))/dR^2;
            e(n,n+1) = Rmesh(1) * 1/dR^2;

            f(n-1) = phi(n-1)/Rmesh(n-1) * (wn-2*w(n-1)+w(n-2))/dR^2 + 1/Rmesh(n-1) * (phin-phi(n-2))/(2*dR) * (wn-w(n-2))/(2*dR) + p(i) + 8*lambda_E^2/3*((w_center(i)+H(j)+1+Rmesh(n-1)^2/2-w(n-1))^(-3)-(w_center(i)+H(j)+1+Rmesh(n-1)^2/2-w(n-1))^(-9));
            g(n-1) = (phin-phi(n-1))/dR - phi(n-1)/Rmesh(n-1) + Rmesh(n-1) * (phin-2*phi(n-1)+phi(n-2))/dR^2 + 1/2 * ( (wn-w(n-1))/dR )^2;

            e(n-1,n-2) = phi(n-1)/Rmesh(n-1) * 1/dR^2 + 1/Rmesh(n-1) * (phin-phi(n-2))/(2*dR) * (-1)/(2*dR);
            e(n-1,n-1) = phi(n-1)/Rmesh(n-1) * (-2)/dR^2 + 8*lambda_E^2/3*(3*(w_center(i)+H(j)+1+Rmesh(n-1)^2/2-w(n-1))^(-4)-9*(w_center(i)+H(j)+1+Rmesh(n-1)^2/2-w(n-1))^(-10));
            e(n-1,2*n-3) = 1/Rmesh(n-1) * (-1)/(2*dR) * (wn-w(n-2))/(2*dR);
            e(n-1,2*n-2) = 1/Rmesh(n-1) * (wn-2*w(n-1)+w(n-2))/dR^2 + 1/Rmesh(n-1) * (1/(1-v*dR/Rf))/(2*dR) * (wn-w(n-2))/(2*dR);

            e(2*n-2,n-2) = 0;
            e(2*n-2,n-1) = w(n-1)/dR^2;
            e(2*n-2,2*n-3) = Rmesh(n-1) * 1/dR^2;
            e(2*n-2,2*n-2) = (1/(1-v*dR/Rf)-1)/dR - 1/Rmesh(n-1) + Rmesh(n-1) * (1/(1-v*dR/Rf)-2)/dR^2;
   
            m(1:n-1) = w(1:n-1);
            m(n:2*n-2) = phi(1:n-1);
   
            b(1:n-1) = f(1:n-1);
            b(n:2*n-2) = g(1:n-1);
   
            e = sparse(e);
            delta = -e\b;
   
            R_temp = delta./m;
            iteration = iteration + 1;
   
            m = m + beta * delta;
   
            w(1:n-1) = m(1:n-1);
            phi(1:n-1) = m(n:2*n-2);

            % If the residual increases, exit the loop
            if iteration < iteration_middle
                R = R_temp;
            elseif norm(R_temp)>norm(R)
                break
            end

        end

        for k=1:n-1
            F(j,i) = F(j,i)+16/3*pi*lambda_E^2*Rmesh(k)*dR*((w_center(i)+H(j)+Rmesh(k)^2/2-w(k)+1)^(-3)-(w_center(i)+H(j)+Rmesh(k)^2/2-w(k)+1)^(-9));
        end
        F(j,i) = F(j,i)+16/3*pi*lambda_E^2*Rf*dR*((w_center(i)+H(j)+Rf^2/2-wn+1)^(-3)-(w_center(i)+H(j)+Rf^2/2-wn+1)^(-9));

        if j>1 && F(j,i)>0 && F(j,i)<F(j-1,i)
            Fmax(i) = F(j-1,i);

            % if Fmax(i)/pi/lambda_E^2 >= 1
                H_pulloff(i) = H(j-1);
                w_pulloff(:,i) = w_temp;
                phi_pulloff(:,i) = phi_temp;

                flag = max(j-500,0);
        
                if (iteration == iteration_max) || ( iteration > iteration_middle && norm(R_temp)>norm(R) )
                    warning_index(i) = 1;
                end

                break
            % end
        end
    
        w_temp = w;
        phi_temp = phi;

    end
end

filename = ['lambda_E_',num2str(lambda_E),'_Rf_',num2str(Rf),'_Rs_',num2str(Rs),'.mat'];
save(filename)

figure_count = 1;

figure(figure_count)
hold on
leg = cell(1,length(lambda));
c = parula(length(lambda));

for i=1:length(p)
    plot(H,F(:,i)/pi/lambda_E^2,'linewidth',1.5,'Color',c(i,:))
    leg{i} = ['$\lambda=$',num2str(lambda(i))];
end

axis on
ax = gca;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.MinorGridAlpha = 0.1;
ax.GridAlpha = 0.1;
ax.TickLength = [0.02 0.02];
%ax.XLim = [-0.3,2];
ax.YLim = [-0.1,2.05];
%ax.XDir = 'reverse';
% legend(leg,'Interpreter','latex','location','best','FontSize',11)
xlabel('$H$','interpreter','latex','fontsize',12)
ylabel('$F/{\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)
title(['$\bar{R_f}=$',num2str(Rf),'$\quad \bar{R_s}=$',num2str(Rs),'$\quad \lambda_E=$',num2str(lambda_E),'$\quad \lambda_T=\infty$'],'Interpreter','latex','fontsize',12)

figure_count = figure_count + 1;

figure(figure_count)
hold on

plot(lambda_p,Fmax/pi/lambda_E^2,'linewidth',1.5,'Color','#E49273')
plot(lambda_p,Fmax/pi/lambda_E^2.*(1+2*w_center/Rf^2),'linewidth',1.5,'Color','#A8D0DB')

xl = xline(lambda_E,'--k','$\lambda_E$','LineWidth',1,'Interpreter','latex','FontSize',12);
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
xl.LabelOrientation = 'horizontal';

line([1e-3,10^(-1.5)],[2,2],'linewidth',2.5,'linestyle','--','color','k')
line([10^(1.5),1e3],[1,1],'linewidth',2.5,'linestyle','--','color','k')

c = [.01 .72 .77
    .99 .49 .00
    .68 .42 .89];

% scatter([1e-2,10^(-0.5),10],[Fmax(11)/pi/lambda(11)^2,Fmax(26)/pi/lambda(26)^2,Fmax(41)/pi/lambda(41)^2],75,c,'Marker','x','LineWidth',2)

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
% ax.YLim = [1,2];
%ax.XDir = 'reverse';
xlabel('$\lambda_p$','interpreter','latex','fontsize',12)
ylabel('$F/{\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)
title(['$\mathcal{R}_f=$',num2str(Rf),'$\quad \mathcal{R}_s=$',num2str(Rs),'$\quad \lambda_E=$',num2str(lambda_E),'$\quad \lambda_T=\infty$'],'Interpreter','latex','fontsize',12)

figure_count = figure_count + 1;

find(warning_index)

toc




