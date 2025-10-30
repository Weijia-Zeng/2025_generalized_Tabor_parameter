clear; clc

global A Rf Rs P Gamma T Adh v

Adh     = 1;        % whether adhesion effects are considered: 1 for yes and 0 for no
P       = 5e-5;     % dimensionless pressure
T       = 0.01;     % dimensionless pretension
Gamma   = 0.05;     % dimensionless adhesion energy
alpha   = 0.05;     % ratio of Rs to Rf
v       = 0.165;    % Poisson ratio

n       = 300;      % number of numerical points for one curve
nmesh   = 1000;     % number of grid nodes (outside the contact area)

temp    = 1/2 * ( sqrt( ( (P/alpha)^(2/3)+T )^2+4*Gamma )-((P/alpha)^(2/3)+T) );    % temporary variable
Rs      = sqrt(1/temp);       % dimensionless sphere radius
Rf      = Rs/alpha;           % dimensionless membrane radius
Rd      = zeros(n,1);         % dimensionless height of the lower point of the sphere
RF      = zeros(n,1);         % dimensionless external force without the consideration of pressure
RF1      = zeros(n,1);        % dimensionless external force with the consideration of pressure
A0      = 0.0001*Rs;          % mimimum dimensionless contact radius (used to approximate cases where the contact radius is zero)
RA      = logspace(log10(A0),log10(Rs/0.9),n);   % range of variation in the dimensionless contact radius


%%%%%%%%%%%%%%%%%%%%%%%%% main code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% calculate bubble height %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        delta0  = -A^2/2+sol.y(1,1);     % dimensionless bubble height
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
        delta0  = -A^2/2+sol.y(1,1);     % dimensionless bubble height
    end

end

ratio1 = delta0/Rf/Rs;          % delta_p/Rf
ratio2 = 2*delta0/Rf^2;         % kappa_of_membrane/(1/Rs)

%%%%%%%%%%%%%%%%% calculate the F-delta and F-A curves %%%%%%%%%%%%%%%

for j = 1:n

    % use the previous results as the initial solution
    A = RA(j);
    xmesh   = logspace(log10(A),log10(Rf),nmesh);
    y       = deval(sol, xmesh);
    solinit.x = xmesh;
    solinit.y = y;

end

sol     = bvp5c(@bvpfcn, @bcfcn, solinit, options);

% if abs(sol.stats.maxerr) > 1e-7
%     solinit = bvpinit(xmesh, @guess);
%     sol     = bvp5c(@bvpfcn, @bcfcn, solinit, options);
% end

sol     = bvp5c(@bvpfcn, @bcfcn, sol, options);
sol     = bvp5c(@bvpfcn, @bcfcn, sol, options);

RF(j)   = -2*sol.y(2,1)*sol.y(3,1)/Gamma/Rs^2;
RF1(j)  = -2*sol.y(2,1)*sol.y(3,1)/Gamma/Rs^2-P*A^2/Gamma/Rs^2;
Rd(j)   = -A^2/2+sol.y(1,1);
j

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

filename = ['JKR_limit_P_',num2str(P),'_T_',num2str(T),'_Gamma_',num2str(Gamma),'_alpha_',num2str(alpha),'.mat'];
save(filename)

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

