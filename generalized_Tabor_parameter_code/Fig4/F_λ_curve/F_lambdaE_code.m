n = 15000;                    % grid number
start = -0.1;                 % starting position of the indenter: (H_start-z0)/z0
cutoff = 1000;                % stop position
H = linspace(start,cutoff,(cutoff-start)*100+1);
v = 0.165;

Rf = 500;                     % dimensionless membrane radius
Rs = 20;                      % dimensionless sphere radius
start = -3;
cutoff = 2;
lambda_E = logspace(start,cutoff,(cutoff-start)*10+1);
lambda = lambda_E;

% initial guess (n-1 components)
W = 0.005*Rf*(1-((1:n-1)'./n).^2);           
PHI = 0.001*ones(n-1,1);

beta = 1;                     % relaxation factor


%%%%%%%%%%%%%%%%  start solving %%%%%%%%%%%%%%%%%%%

tic
f = zeros(n-1,1);
g = zeros(n-1,1);
e = zeros(2*n-2,2*n-2); % Jacobi matrix
b = zeros(2*n-2,1);     % formed by connecting vectors f and g
m = zeros(2*n-2,1);     % formed by connecting vectors w and phi

% record the conditions of pull-off points
w_pulloff = zeros(n-1,length(lambda));
phi_pulloff = zeros(n-1,length(lambda));
H_pulloff = zeros(1,length(lambda));

dR = Rf/n;
Rmesh = (1:n-1).*dR;
     
F = zeros(length(H),length(lambda));
Fmax = zeros(1,length(lambda));


%%%%%%%%%%%%%%%%%%%% Lifting the sphere %%%%%%%%%%%%%%%%%%%%%%%

flag = 0;
iteration_max = 1000;
iteration_middle = iteration_max/2;
warning_index = zeros(1,length(lambda_E));  % Indicator for curve truncation due to iteration termination

for i = 1:length(lambda_E)

    if i == 1
        w = W;
        phi = PHI;
    end
  
    for j=1:length(H)
      
        R = ones(2*n-2,1);        % relative incremental step
        iteration = 0;
        if j < flag
            continue
        end
    
        % start to solve
        while(norm(R)>1e-8 && iteration ~= iteration_max)
    
            % avoid the membrane penetrateing the sphere surface
            for k=1:n-1
                if w(k) >= H(j)+0.98+1/2*Rmesh(k)^2
                    w(k) = H(j)+0.98+1/2*Rmesh(k)^2;
                end
            end

            % boundary conditions
            w0 = (4*w(1)-w(2))/3;
            wn = 0;
            phi0 = phi(2)-2*v*dR*phi(1)/Rmesh(1);
            phin = ( 4*phi(n-1)-phi(n-2) ) / ( 3-2*v*dR/Rf );

            for k=2:(n-2) 
                f(k) = phi(k)/Rmesh(k) * (w(k+1)-2*w(k)+w(k-1))/dR^2 + 1/Rmesh(k) * (phi(k+1)-phi(k-1))/(2*dR) * (w(k+1)-w(k-1))/(2*dR) + 8*lambda_E(i)^2/3*((H(j)+1+Rmesh(k)^2/2-w(k))^(-3)-(H(j)+1+Rmesh(k)^2/2-w(k))^(-9));
                g(k) = (phi(k+1)-phi(k-1))/(2*dR) - phi(k)/Rmesh(k) + Rmesh(k) * (phi(k+1)-2*phi(k)+phi(k-1))/dR^2 + 1/2 * ( (w(k+1)-w(k-1))/(2*dR) )^2;

                e(k,k-1) = phi(k)/Rmesh(k) * 1/dR^2 - 1/Rmesh(k) * (phi(k+1)-phi(k-1))/(2*dR) * 1/(2*dR);
                e(k,k) = phi(k)/Rmesh(k) * (-2)/dR^2 + 8*lambda_E(i)^2/3*(3*(H(j)+1+Rmesh(k)^2/2-w(k))^(-4)-9*(H(j)+1+Rmesh(k)^2/2-w(k))^(-10));
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

            f(1) = phi(1)/Rmesh(1) * (w(2)-2*w(1)+w0)/dR^2 + 1/Rmesh(1) * (phi(2)-phi0)/(2*dR) * (w(2)-w0)/(2*dR) + 8*lambda_E(i)^2/3*((H(j)+1+Rmesh(1)^2/2-w(1))^(-3)-(H(j)+1+Rmesh(1)^2/2-w(1))^(-9));
            g(1) = (phi(2)-phi0)/(2*dR) - phi(1)/Rmesh(1) + Rmesh(1) * (phi(2)-2*phi(1)+phi0)/dR^2 + 1/2 * ( (w(2)-w0)/(2*dR) )^2;

            e(1,1) =  phi(1)/Rmesh(1) * (-2+4/3)/dR^2 + 1/Rmesh(1) * (phi(2)-phi0)/(2*dR) * (-4/3)/(2*dR) + 8*lambda_E(i)^2/3*(3*(H(j)+1+Rmesh(1)^2/2-w(1))^(-4)-9*(H(j)+1+Rmesh(1)^2/2-w(1))^(-10));
            e(1,2) = phi(1)/Rmesh(1) * (1-1/3)/dR^2 + 1/Rmesh(1) * (phi(2)-phi0)/(2*dR) * (1+1/3)/(2*dR);
            e(1,n) = 1/Rmesh(1) * (w(2)-2*w(1)+w0)/dR^2 + 1/Rmesh(1) * (2*v*dR/Rmesh(1))/(2*dR) * (w(2)-w0)/(2*dR);
            e(1,n+1) = 1/Rmesh(1) * (1+1)/(2*dR) * (w(2)-w0)/(2*dR);

            e(n,1) = (w(2)-w0)/(2*dR) * (-4/3)/(2*dR);
            e(n,2) = (w(2)-w0)/(2*dR) * (1+1/3)/(2*dR);
            e(n,n) = (2*v*dR/Rmesh(1))/(2*dR) - 1/Rmesh(1) + Rmesh(1) * (-2-2*v*dR/Rmesh(1))/dR^2;
            e(n,n+1) = Rmesh(1) * (1+1)/dR^2;

            f(n-1) = phi(n-1)/Rmesh(n-1) * (wn-2*w(n-1)+w(n-2))/dR^2 + 1/Rmesh(n-1) * (phin-phi(n-2))/(2*dR) * (wn-w(n-2))/(2*dR) + 8*lambda_E(i)^2/3*((H(j)+1+Rmesh(n-1)^2/2-w(n-1))^(-3)-(H(j)+1+Rmesh(n-1)^2/2-w(n-1))^(-9));
            g(n-1) = (phin-phi(n-2))/(2*dR) - phi(n-1)/Rmesh(n-1) + Rmesh(n-1) * (phin-2*phi(n-1)+phi(n-2))/dR^2 + 1/2 * ( (wn-w(n-2))/(2*dR) )^2;

            e(n-1,n-2) = phi(n-1)/Rmesh(n-1) * 1/dR^2 + 1/Rmesh(n-1) * (phin-phi(n-2))/(2*dR) * (-1)/(2*dR);
            e(n-1,n-1) = phi(n-1)/Rmesh(n-1) * (-2)/dR^2 + 8*lambda_E(i)^2/3*(3*(H(j)+1+Rmesh(n-1)^2/2-w(n-1))^(-4)-9*(H(j)+1+Rmesh(n-1)^2/2-w(n-1))^(-10));
            e(n-1,2*n-3) = 1/Rmesh(n-1) * (-1/(3-2*v*dR/Rf)-1)/(2*dR) * (wn-w(n-2))/(2*dR);
            e(n-1,2*n-2) = 1/Rmesh(n-1) * (wn-2*w(n-1)+w(n-2))/dR^2 + 1/Rmesh(n-1) * (4/(3-2*v*dR/Rf))/(2*dR) * (wn-w(n-2))/(2*dR);

            e(2*n-2,n-2) = w(n-2)/(2*dR)^2;
            e(2*n-2,n-1) = 0;
            e(2*n-2,2*n-3) = (-1/(3-2*v*dR/Rf)-1)/(2*dR) + Rmesh(n-1) * (-1/(3-2*v*dR/Rf)+1)/dR^2;
            e(2*n-2,2*n-2) = (4/(3-2*v*dR/Rf))/(2*dR) - 1/Rmesh(n-1) + Rmesh(n-1) * (4/(3-2*v*dR/Rf)-2)/dR^2;
   
            m(1:n-1) = w(1:n-1);
            m(n:2*n-2) = phi(1:n-1);
   
            b(1:n-1) = f(1:n-1);
            b(n:2*n-2) = g(1:n-1);
   
            e = 100000*e;
            e = sparse(e);
            delta = -e\(100000*b);
   
            R_temp = delta./m;    % norm of the relative incremental step
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
            F(j,i) = F(j,i)+16/3*pi*lambda_E(i)^2*Rmesh(k)*dR*((H(j)+Rmesh(k)^2/2-w(k)+1)^(-3)-(H(j)+Rmesh(k)^2/2-w(k)+1)^(-9));
        end
        F(j,i) = F(j,i)+16/3*pi*lambda_E(i)^2*Rf*dR*((H(j)+Rf^2/2-wn+1)^(-3)-(H(j)+Rf^2/2-wn+1)^(-9));

        if j>1 && F(j,i)>0 && F(j,i)<F(j-1,i)
            Fmax(i) = F(j-1,i);

            if Fmax(i)/pi/lambda_E(i)^2 >= 1
                H_pulloff(i) = H(j-1);
                w_pulloff(:,i) = w_temp;
                phi_pulloff(:,i) = phi_temp;

                flag = j-100;
        
                if (iteration == iteration_max) || ( iteration > iteration_middle && norm(R_temp)>norm(R) )
                    warning_index(i) = 1;
                end

                break
            end
        end
    
        w_temp = w;
        phi_temp = phi;

    end
end

filename = 'lambda_E_data.mat';
save(filename)

figure_count = 1;

%% F-delta curve
figure(figure_count)
hold on
leg = cell(1,length(lambda));
c = parula(length(lambda));

for i=1:length(lambda)
    plot(H,F(:,i)/pi/lambda_E(i)^2,'linewidth',1.5)
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
legend(leg,'Interpreter','latex','location','best','FontSize',11)
xlabel('$H$','interpreter','latex','fontsize',12)
ylabel('$F/{\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)
title(['$\bar{R_f}=$',num2str(Rf),'$\quad \bar{R_s}=$',num2str(Rs),'$\quad \lambda_p=\infty \quad \lambda_T=\infty$'],'Interpreter','latex','fontsize',12)

figure_count = figure_count + 1;

%% Fmax-lambda curve
figure(figure_count)
hold on

plot(lambda_E,Fmax/pi./lambda_E.^2,'linewidth',2,'Color','#E49273')

line([1e-3,10^(-1.5)],[2,2],'linewidth',2.5,'linestyle','--','color','k')
line([10^(0.5),1e2],[1,1],'linewidth',2.5,'linestyle','--','color','k')

c = [.01 .72 .77
    .99 .49 .00
    .68 .42 .89];

scatter([1e-2,10^(-0.5),10],[Fmax(11)/pi/lambda(11)^2,Fmax(26)/pi/lambda(26)^2,Fmax(41)/pi/lambda(41)^2],75,c,'Marker','x','LineWidth',2)

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
title(['$\bar{R_f}=$',num2str(Rf),'$\quad \bar{R_s}=$',num2str(Rs),'$\quad \lambda_p=\infty \quad \lambda_T=\infty$'],'Interpreter','latex','fontsize',12)

figure_count = figure_count + 1;

find(warning_index)

toc




