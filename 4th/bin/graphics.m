folder1 = ['results/'];
folders = ['graphics/'];
        
color1=['#c888a9'];
color2=['#a9c888'];
color3=['#88a9c8'];

X = load([folder1 'x_i.txt']);
XX = load([folder1 'xx_i.txt']);
T = load([folder1 'T.txt']);
Lambda = load([folder1 'lambda_i.txt']);
y=zeros([1 35]);

v=figure;

plot(X,T,'-*',...
    'LineWidth',2,...
    'Color', color1,...
    'MarkerSize',2,... 
    'MarkerFaceColor',color1);
title('T(x)', 'fontsize', 12);
xlabel('x, м', 'fontsize', 12);
ylabel('T, K', 'fontsize', 12);
xlim([-X(3),15e-3]);
saveas(v, [folders 'T_x.png']);
close


f=figure;

plot(X,Lambda,'-*',...
    'LineWidth',2,...
    'Color', color1,...
    'MarkerSize',2,... 
    'MarkerFaceColor',color1);
title('Lambda(x)', 'fontsize', 12);
xlabel('x, м', 'fontsize', 12);
ylabel('Lambda, Вт/(м*K)', 'fontsize', 12);
xlim([-X(3),15e-3]);
saveas(f, [folders 'lambda_x.png']);
close


g=figure;

plot(X,y,'*',...
    'LineWidth',2,...
    'Color', color1,...
    'MarkerSize',2,... 
    'MarkerFaceColor',color1);
hold on;
plot(XX,y,'|',...
    'LineWidth',2,...
    'Color', color2,...
    'MarkerSize',6,... 
    'MarkerFaceColor',color2);
hold off;
legend('x_i','xx_i')
xlabel('x, м', 'fontsize', 12);
xlim([0.0132,0.01504]);
set(gca,'ytick',[])
saveas(g, [folders 'x_xx.png']);
close