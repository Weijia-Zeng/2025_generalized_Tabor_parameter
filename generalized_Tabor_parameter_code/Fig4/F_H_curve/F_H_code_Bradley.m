n = 6000;
start = -0.2;
cutoff = 650;
H = linspace(start,cutoff,(cutoff-start)*100+1);
v = 0.165;

Rf = 100;
Rs = 20;
lambda_E = [1e-3 1e-2 1e-1];
lambda = lambda_E;

W = 0.005*Rf*(1-((1:n-1)'./n).^2);            
PHI = 0.001*ones(n-1,1);
beta = 1;


%%%%%%%%%%%%%%%%  start solving  %%%%%%%%%%%%%%%%%%%

tic
f = zeros(n-1,1);
g = zeros(n-1,1);
e = zeros(2*n-2,2*n-2);
b = zeros(2*n-2,1);
m = zeros(2*n-2,1);

w_pulloff = zeros(n-1,length(lambda));
phi_pulloff = zeros(n-1,length(lambda));
H_pulloff = zeros(1,length(lambda));

w_initial = zeros(n-1,length(lambda));
phi_initial = zeros(n-1,length(lambda));

dR = Rf/n;
Rmesh = (1:n-1).*dR;
     
F_1 = zeros(length(H),length(lambda));
F_2 = zeros(length(H),length(lambda));
Fmax = zeros(1,length(lambda));


%%%%%%%%%%%%%%%%%%%%% approaching the sphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iteration_max = 1000;
iteration_middle = iteration_max;

for i = 1:length(lambda)

    w = W;
    phi = PHI;
  
    for j = length(H):-1:1
      
        R = ones(2*n-2,1);
        iteration = 0;
    
        % start to solve
        while(norm(R)>1e-8 && iteration ~= iteration_max)

            % avoid the membrane penetrateing the sphere surface
            for k=1:n-1
                if w(k)>=H(j)+0.98+1/2*Rmesh(k)^2
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
            F_1(j,i) = F_1(j,i)+16/3*pi*lambda_E(i)^2*Rmesh(k)*dR*((H(j)+Rmesh(k)^2/2-w(k)+1)^(-3)-(H(j)+Rmesh(k)^2/2-w(k)+1)^(-9));
        end
        F_1(j,i) = F_1(j,i)+16/3*pi*lambda_E(i)^2*Rf*dR*((H(j)+Rf^2/2-wn+1)^(-3)-(H(j)+Rf^2/2-wn+1)^(-9));

    end

    w_initial(:,i) = w;
    phi_initial(:,i) = phi;
end


%%%%%%%%%%%%%%%%%%% Lifting the sphere %%%%%%%%%%%%%%%%%%

for i=1:length(lambda_E)

    w = w_initial(:,i);
    phi = phi_initial(:,i);

    mark = 0;      % refresh variable 'mark' to indicate whether the pull-off force has been recorded
  
    for j=1:length(H)
      
        R = ones(2*n-2,1);
        iteration = 0;
    
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
   
            R_temp = delta./m;
            iteration = iteration + 1;
   
            m = m + beta*delta;
   
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
            F_2(j,i) = F_2(j,i)+16/3*pi*lambda_E(i)^2*Rmesh(k)*dR*((H(j)+Rmesh(k)^2/2-w(k)+1)^(-3)-(H(j)+Rmesh(k)^2/2-w(k)+1)^(-9));
        end
        F_2(j,i) = F_2(j,i)+16/3*pi*lambda_E(i)^2*Rf*dR*((H(j)+Rf^2/2-wn+1)^(-3)-(H(j)+Rf^2/2-wn+1)^(-9));

        if j>1 && F_2(j,i)<F_2(j-1,i) && mark == 0
            Fmax(i) = F_2(j-1,i);
            mark = 1;

            H_pulloff(i) = H(j-1);
            w_pulloff(:,i) = w_temp;
            phi_pulloff(:,i) = phi_temp;
        end

        w_temp = w;
        phi_temp = phi;

    end
end


%%%%%%%%%%%%%%%%%%% Bradley curve %%%%%%%%%%%%%%%%%%%%%%

F_Bradley = zeros(length(H),1);

for i=1:length(H)

    for k=1:n-1
        F_Bradley(i) = F_Bradley(i)+16/3*Rmesh(k)*dR*((H(i)+Rmesh(k)^2/2+1)^(-3)-(H(i)+Rmesh(k)^2/2+1)^(-9));
    end
    F_Bradley(i) = F_Bradley(i)+16/3*Rf*dR*((H(i)+Rf^2/2+1)^(-3)-(H(i)+Rf^2/2+1)^(-9));

end

%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%

figure(1)
hold on
leg = {'$Bradley\ limit$','$\lambda_E=10^{-3}$','$\lambda_E=10^{-2}$','$\lambda_E=10^{-1}$'};
c = [.01 .72 .77
    .99 .49 .00
    .68 .42 .89];
Index = zeros(2,length(lambda));
one_line = zeros(1,length(lambda));

plot(H,F_Bradley,'color','k','LineWidth',2.5)

for i=1:length(lambda)

    A = isfinite(F_1(:,i));         % Check 'NaN' and 'inf' components
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
        plot(H,F_1(:,i)/(pi*lambda_E(i)^2),'color',c(i,:),'LineWidth',2)
    else
        plot(H(Index(1,i)+1:end),F_1(Index(1,i)+1:end,i)/(pi*lambda_E(i)^2),'color',c(i,:),'linewidth',2)
    end
end

for i=1:length(lambda)
    if one_line(i) == 1
        continue
    else
        plot(H(1:Index(2,i)),F_2(1:Index(2,i),i)/pi/lambda_E(i)^2,'color',c(i,:),'linewidth',2)
    end
end

for i=1:length(lambda)
    if one_line(i) == 1
        continue
    else
        line([H(Index(1,i)),H(Index(1,i)+1)],[F_2(Index(1,i),i)/pi/lambda_E(i)^2,F_1(Index(1,i)+1,i)/pi/lambda_E(i)^2],'color',c(i,:),'linestyle','--')
        line([H(Index(2,i)),H(Index(2,i)+1)],[F_2(Index(2,i),i)/pi/lambda_E(i)^2,F_1(Index(2,i)+1,i)/pi/lambda_E(i)^2],'color',c(i,:),'linestyle','--')
    end
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
%ax.XLim = [-0.3,1.8];
ax.YLim = [-0.1,2.05];
%ax.XDir = 'reverse';
legend(leg,'Interpreter','latex','location','best','FontSize',12.5)
xlabel('$H$','interpreter','latex','fontsize',12)
ylabel('$F/{2\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)
title(['$\bar{R_f}=$',num2str(Rf),'$\quad \bar{R_s}=$',num2str(Rs),'$\quad \lambda_p=\infty \quad \lambda_T=\infty$'],'Interpreter','latex','fontsize',12)

filename = 'F_H_data_Bradley.mat';
% filename = 'F_H_data_approach_JKR.mat';
save(filename)

toc


