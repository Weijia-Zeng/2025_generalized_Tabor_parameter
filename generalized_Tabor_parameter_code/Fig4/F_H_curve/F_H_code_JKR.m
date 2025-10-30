clear; clc

global A Rf Rs P Gamma T Adh v

% % read numerical solutions as initial guess for new lamdba_E
% filename = 'F_H_data_JKR_lambdaE_100.mat';
% load(filename)

% microscopic dimensionless system values
lambda_E = 10;
Rf_micro = 100;
Rs_micro = 20;

% macroscopic dimensionless system values
Adh     = 1;
P       = 0;
T       = 0;
Gamma   = lambda_E^2/Rs_micro^4;
alpha   = Rs_micro/Rf_micro;
v       = 0.165;
delta0  = 0;

n       = 300;
nmesh   = 1000;

temp    = 1/2 * ( sqrt( ( (P/alpha)^(2/3)+T )^2+4*Gamma )-((P/alpha)^(2/3)+T) );
Rs      = sqrt(1/temp);
Rf      = Rs/alpha;
Rd      = zeros(n,1);
RF      = zeros(n,1);
RF1      = zeros(n,1);
A0      = 0.00015*Rs;
RA      = logspace(log10(A0),log10(Rs/2.7),n);

initial_switch = 0;   % whther to use numerical solutions as the inital guess (1 for yes and 0 for no) (parameter approximation method)

if initial_switch ~= 1
    % record the initial solution for handling computationally challenges
    w_initial = zeros(nmesh,n);
    dw_initial = zeros(nmesh,n);
    phi_initial = zeros(nmesh,n);
    dphi_initial = zeros(nmesh,n);
end

%%%%%%%%%%%%%%%%%%%%%%%%% main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% calculate bubble %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if P ~= 0

    if Adh    == 1

        Adh     = 0;
        A       = A0;
        xmesh   = logspace(log10(A),log10(Rf),nmesh);
        solinit = bvpinit(xmesh, @guess);
        options = bvpset('RelTol',1e-9,'Stats','on');
        sol     = bvp5c(@bvpfcn, @bcfcn, solinit, options);
        sol     = bvp5c(@bvpfcn, @bcfcn, sol, options);
        sol     = bvp5c(@bvpfcn, @bcfcn, sol, options);
        delta0  = -A^2/2+sol.y(1,1);
        Adh     = 1;
    end

    if Adh    == 0
        A       = A0;
        xmesh   = logspace(log10(A),log10(Rf),nmesh);
        solinit = bvpinit(xmesh, @guess);
        options = bvpset('RelTol',1e-9,'Stats','on');
        sol     = bvp5c(@bvpfcn, @bcfcn, solinit, options);
        sol     = bvp5c(@bvpfcn, @bcfcn, sol, options);
        sol     = bvp5c(@bvpfcn, @bcfcn, sol, options);
        delta0  = -A^2/2+sol.y(1,1);
    end

    ratio1 = delta0/Rf/Rs;
    ratio2 = 2*delta0/Rf^2;

end

%%%%%%%%%%%%%%%%%%%% start to solve %%%%%%%%%%%%%%%%%%%%%%

options = bvpset('RelTol',1e-9,'Stats','on');

for j = 1:n

    A = RA(j);
    xmesh   = logspace(log10(A),log10(Rf),nmesh);

    if initial_switch == 1
        solinit.x = xmesh;
        solinit.y = [w_initial(:,j)'; dw_initial(:,j)'; phi_initial(:,j)'; dphi_initial(:,j)'];
    else
        if j ~= 1 
            y       = deval(sol, xmesh);
            solinit.x = xmesh;
            solinit.y = y;
        else
            solinit = bvpinit(xmesh,@guess);
        end
    end

    sol     = bvp5c(@bvpfcn, @bcfcn, solinit, options);

    if abs(sol.stats.maxerr) > 1e-5
        solinit = bvpinit(xmesh, @guess);
        sol     = bvp5c(@bvpfcn, @bcfcn, solinit, options);
    end

    sol     = bvp5c(@bvpfcn, @bcfcn, sol, options);
    sol     = bvp5c(@bvpfcn, @bcfcn, sol, options);

    RF(j)   = -2*sol.y(2,1)*sol.y(3,1)/Gamma/Rs^2;
    RF1(j)  = -2*sol.y(2,1)*sol.y(3,1)/Gamma/Rs^2-P*A^2/Gamma/Rs^2;
    Rd(j)   = -A^2/2+sol.y(1,1);
    j

    % record the numerical solutions for initial guess
    y = deval(sol, xmesh);
    w_initial(:,j) = y(1,:);
    dw_initial(:,j) = y(2,:);
    phi_initial(:,j) = y(3,:);
    dphi_initial(:,j) = y(4,:);
end

figure(1)
hold on

% plot((Rd-delta0),RF,'LineWidth',1.5)
plot((Rd-delta0),RF1,'LineWidth',1.5)

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
% ax.XLim = [0,5];
% ax.YLim = [0,2];
%ax.XDir = 'reverse';
xlabel('$\Delta=\delta/\delta_*$','interpreter','latex','fontsize',12)
ylabel('$F/{\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)

figure(2)
hold on

plot(RA,RF1,'LineWidth',1.5)

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
% ax.XLim = [0,5];
% ax.YLim = [0,2];
%ax.XDir = 'reverse';
legend('$\rm JKR$','Interpreter','latex','location','best','FontSize',11)
xlabel('$A=a/a_*$','interpreter','latex','fontsize',12)
ylabel('$F/{\pi\gamma{R_s}}$','interpreter','latex','fontsize',12)

filename = ['F_H_data_JKR_lambdaE_',num2str(lambda_E),'.mat'];
% save(filename)

function g = guess(x)
global Rf T

    g = [0.01*(1-x.^2/Rf^2)
        -2*0.01*x/Rf^2
        (T+0.001).*x
        (T+0.001)];
    
end

function dydx = bvpfcn(x,y)
global P Rs

%y(1)=w;y(2)=w';y(3)=phi;y(4)=phi'

if P ~=0
    dydx = [y(2)
            -x*P/y(3)-y(2)*y(4)/y(3)
            y(4)
            -y(4)/x+y(3)/x/x-1/2/x*y(2)^2/Rs^2];
else
    dydx = [y(2)
            -y(2)*y(4)/y(3)
            y(4)
            -y(4)/x+y(3)/x/x-1/2/x*y(2)^2/Rs^2];
end
end

function res = bcfcn(ya,yb) % boundary conditions
global A Gamma Adh v T Rf Rs
res = [A-ya(2)-(2*A*Gamma/ya(3))^0.5*Rs*Adh
       ya(3)-A*ya(4)-A^3/8/Rs^2
       yb(1)
       yb(4)-v*yb(3)/Rf-(1-v)*T];
end

