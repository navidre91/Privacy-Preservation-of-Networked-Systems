clear
clc

colors=[1 0.5 0;1 0 1;0 1 1;1 0 0;0 1 0;0 0 1; 0 0 0;0.25 1 0.25];
colors=[[0,0,0]+0.10;[0,0,0]+0.20;[0,0,0]+0.30;[0,0,0]+0.40;[0,0,0]+0.50;[0,0,0]+0.60; [0,0,0]+0.70;[0,0,0]+0.80];
colors=[[0,0,0];[0,0,0]+0.20;[0,0,0]+0.40;[0,0,0]+0.60;[0,0,0]+0.80];

x=[3;2;5;-3;-1];

dt=0.00001;

dout=[3;1;1;2;2];

T=15;


A=[0 1 0 1 1;0 0 1 0 0;1 0 0 0 0;1 0 0 0 1;1 0 0 1 0];
L=diag(dout)-A;
i=1;

x_hat_ext=0;
x_hat_2_ext=0;
x_hat_4=0;
x_hat_5=0;
z_hat_2_ext=0;
z_hat_4=0;
z_hat_5=0;

g_rec=[];
for t=0:dt:T
    
    g=[];
    f=[];
    
    for j=1:5
    g=[g;sin(j*pi/12+j*pi*t^2)];
    f=[f;-dout(j)*(sin(j*pi/12)+cos(j*pi/12))*sqrt(2*j)/4/j*exp(-t)];
    end
        
    y=x+g;
    z_hat_2_ext=x_hat_2_ext+x_hat_ext;
    z_hat_4=x_hat_4+x(1);
    z_hat_5=x_hat_5+x(1);
    
    z_hat_2_ext_rec(i)=z_hat_2_ext;
    z_hat_4_rec(i)=z_hat_4;
    z_hat_5_rec(i)=z_hat_5;
    y_rec1(:,i)=y;
    x_rec1(:,i)=x;
    
    i=i+1;
    x_dot=-L*y+f+diag(dout)*g;
    
    x_hat_2_ext_dot=y(2)-y(3);
    x_hat_ext_dot=-x_hat_ext+y(2);
    x_hat_4_dot=(y(4)-y(1))+(y(4)-y(5));
    x_hat_5_dot=(y(5)-y(1))+(y(5)-y(4));
    
    x_hat_2_ext=x_hat_2_ext+x_hat_2_ext_dot*dt;
    x_hat_ext=x_hat_ext+x_hat_ext_dot*dt;
    x_hat_4=x_hat_4+x_hat_4_dot*dt;
    x_hat_5=x_hat_5+x_hat_5_dot*dt;
    x=x+x_dot*dt;
    
    
    
end

tt=0:dt:T;




x=[3;1;6;-3;-1];


A=[0 1 0 1 1;0 0 1 0 0;1 0 0 0 0;1 0 0 0 1;1 0 0 1 0];
L=diag(dout)-A;
i=1;

for t=0:dt:T
    
    g=[];
    f=[];
    
    for j=1:5
    g=[g;sin(j*pi/12+j*pi*t^2)];
    f=[f;-dout(j)*(sin(j*pi/12)+cos(j*pi/12))*sqrt(2*j)/4/j*exp(-t)];
    end
    
    g(2)=g(2)-exp(-1*t)*(-1);
    f(2)=f(2)-1*exp(-1*t)*1;
    
    y=x+g;
    
    y_rec2(:,i)=y;
    
    x_rec2(:,i)=x;
    i=i+1;
    x_dot=-L*y+f+diag(dout)*g;
    x=x+x_dot*dt;
end

t=0:dt:T;

figure('Units','inches',...
    'Position',[5 5 2.5 2],...
    'PaperPositionMode','auto');

for j=1:5
    hold on
    plot(t,x_rec1(j,:),'-','LineWidth',0.4,'Color',colors(j,:))
    plot(t,x_rec2(j,:),'--','LineWidth',1.2,'Color',colors(j,:))
end



legend({'$Actual$','$Alternative$'},'Interpreter','latex');
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('$x^i(t)$','interpreter','latex','FontSize',14)
ylim([-3,6])
set(gca,'FontSize',14)
set(gca, 'YGrid', 'on')

set(gca,...
    'Units','normalized',...
    'XTick',0:5:T(end),...
    'Position',[.16 .25 .75 .68],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')

print('G:\My Drive\Privacy Preservation\MATLAB\Journal Figs\1','-depsc2')


figure('Units','inches',...
    'Position',[5 5 2.5 2],...
    'PaperPositionMode','auto');

for j=1:5
    hold on
    if j==1
    plot(t,y_rec1(j,:)-y_rec2(j,:),'-','LineWidth',1,'Color',colors(5,:))
    end
    if j==2
    plot(t,y_rec1(j,:)-y_rec2(j,:),'--','LineWidth',3,'Color',colors(4,:))
    end
    if j==3
    plot(t(1:10000:end),y_rec1(j,(1:10000:end))-y_rec2(j,(1:10000:end)),'.','LineWidth',0.4,'Color',colors(1,:))
    end
    if j==4
    plot(t(1:100000:end),y_rec1(j,(1:100000:end))-y_rec2(j,(1:100000:end)),'-s','LineWidth',0.4,'Color',colors(1,:),'MarkerSize',4)
    end
    if j==5
    plot(t(1+5000:100000:end),y_rec1(j,(1+5000:100000:end))-y_rec2(j,(1+5000:100000:end)),'-d','LineWidth',0.4,'Color',colors(1,:),'MarkerSize',4)
    end
end
legendflex({'$\delta {y^1}$','$\delta {y^2}$','$\delta {y^3}$','$\delta {y^4}$','$\delta {y^5}$'},'Interpreter','latex','ncol',2)
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('$\delta {y^i}(t)$','interpreter','latex','FontSize',14)
set(gca,'FontSize',14)
set(gca, 'YGrid', 'on')


set(gca,...
    'Units','normalized',...
    'XTick',0:5:T(end),...
    'Position',[.16 .25 .75 .68],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')

print('G:\My Drive\Privacy Preservation\MATLAB\Journal Figs\2','-depsc2')

figure('Units','inches',...
    'Position',[5 5 2.5 2],...
    'PaperPositionMode','auto');
hold on
plot(t,z_hat_4_rec,'--','LineWidth',0.4,'Color',colors(1,:))
plot(t,z_hat_5_rec,'LineWidth',0.4,'Color',colors(4,:))

legend({'$\nu^4$','$\nu^5$'},'interpreter','latex')
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('$\nu(t)$','interpreter','latex','FontSize',14)
set(gca,'FontSize',14)
set(gca, 'YGrid', 'on')


set(gca,...
    'Units','normalized',...
    'YTick',-3:2:5,...
    'XTick',0:5:T(end),...
    'Position',[.16 .25 .75 .68],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')

print('G:\My Drive\Privacy Preservation\MATLAB\Journal Figs\3','-depsc2')

figure('Units','inches',...
    'Position',[5 5 2.5 2],...
    'PaperPositionMode','auto');
hold on
plot(t,z_hat_2_ext_rec,'LineWidth',0.4,'Color',colors(1,:))

% legend({'$\hat{z}^2$'},'interpreter','latex')
xlabel('$t$','interpreter','latex','FontSize',14)
ylabel('$\nu^2(t)$','interpreter','latex','FontSize',14)
set(gca,'FontSize',14)
set(gca, 'YGrid', 'on')


set(gca,...
    'Units','normalized',...
    'XTick',0:5:T(end),...
    'Position',[.16 .25 .75 .68],...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',12,...
    'FontName','Times')

print('G:\My Drive\Privacy Preservation\MATLAB\Journal Figs\4','-depsc2')