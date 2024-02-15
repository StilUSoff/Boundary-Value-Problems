folder1 = ['res/coord/'];
folder2 = ['res/time/'];
folders = ['graphics/'];
        
color1=['#c888a9'];
color2=['#a9c888'];
color3=['#88a9c8'];

C=load([folder1 'bi_alpha.txt']);
AlphaCoord=C(:,1);
BiCoord=C(:,2);
D=load([folder2 'bi_alpha.txt']);
AlphaTime=D(:,1);
BiTime=D(:,2);



n=length(AlphaCoord);


for j=1:1:n
v=figure;
lop=round(AlphaCoord(j));
pup=BiCoord(j);
for i=1:1:3

    A = load([folder1 'alpha' int2str(lop) '_point' int2str(i) '_temp.txt']);
    B = load([folder1 'alpha' int2str(lop) '_point' int2str(i) '_coord.txt']);

    x = B(:,1);
    y = A(:,1);
    if i==1
        plot(x,y,'-*',...
            'LineWidth',1,...
            'Color', color1,...
            'MarkerSize',1,... 
            'MarkerFaceColor',color1);
    elseif i==2
        plot(x,y,'-*',...
            'LineWidth',1,...
            'Color', color2,...
            'MarkerSize',1,... 
            'MarkerFaceColor',color2);
    else
        plot(x,y,'-*',...
            'LineWidth',1,...
            'Color', color3,...
            'MarkerSize',1,... 
            'MarkerFaceColor',color3);
    end
hold on;
end
%title(['T(x), α=' int2str(lop) ', Bi=' num2str(pup)], 'fontsize', 12);
title(['Bi=' num2str(pup)], 'fontsize', 12);
xlabel('x, м', 'fontsize', 12);
ylabel('T, K', 'fontsize', 12);
xlim([-0.01*x(100),x(100)+0.01*x(100)]);
ylim([400,800]);
legend('t=0.1∙Fo_m_a_x', 't=0.5∙Fo_m_a_x', 't=0.9∙Fo_m_a_x','Location','northeast');
hold off;
saveas(v, [folders 'temp_x_alpha' int2str(lop) '.png']);
close
end

n=length(AlphaTime);

for j=1:1:n
k=figure;
lop=round(AlphaTime(j));
pup=BiTime(j);
FP = load([folder2 'fo_points_alpha' int2str(lop) '.txt']);
fnum = FP(:,1);
fpoints = FP(:,2);
fzero = FP(:,3);
for i=1:1:3

C = load([folder2 'alpha' int2str(lop) '_point' int2str(i) '_temp.txt']);
D = load([folder2 'alpha' int2str(lop) '_point' int2str(i) '_time.txt']);

x = D(:,1);
y = C(:,1);
if i==1
        plot(x,y,'-*',...
            'LineWidth',1,...
            'Color', color1,...
            'MarkerSize',1,... 
            'MarkerFaceColor',color1);
elseif i==2
        plot(x,y,'-*',...
            'LineWidth',1,...
            'Color', color2,...
            'MarkerSize',1,... 
            'MarkerFaceColor',color2);
else
        plot(x,y,'-*',...
            'LineWidth',1,...
            'Color', color3,...
            'MarkerSize',1,... 
            'MarkerFaceColor',color3);
end
hold on;
end
plot(fpoints,fzero,'|', ... 
    'MarkerSize',10);
hold on;
%title(['T(t), α=' int2str(lop) ', Bi=' num2str(pup)], 'fontsize', 12);
title(['Bi=' num2str(pup)], 'fontsize', 12);
xlabel('t, c', 'fontsize', 12);
ylabel('T, K', 'fontsize', 12);
xlim([-0.01*x(100),x(100)+0.01*x(100)]);
ylim([400,800]);
legend('x=0', 'x=d/4', 'x=d/2','Location','northeast');
hold off;
saveas(k, [folders 'temp_t_alpha' int2str(lop) '.png']);
close
end