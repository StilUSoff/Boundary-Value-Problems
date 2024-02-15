folder1 = ['res/Im=101/'];
folders = ['graphics/'];
        
color1=['#c888a9'];
color2=['#a9c888'];
color3=['#88a9c8'];

C=load([folder1 'bi_iFORcoord_iFORtime_alpha.txt']);
AlphaCoord=C(:,4);
BiCoord=C(:,1);
x_lim=C(:,2);
time_lim=C(:,3);



n=length(AlphaCoord);


for j=1:1:n
v=figure;
lop=round(AlphaCoord(j));
pup=BiCoord(j);
for i=1:1:3

    A = load([folder1 'alpha' int2str(lop) '_temp_x.txt']);
    B = load([folder1 'alpha' int2str(lop) '_coord.txt']);

    x = B(:,1);
    y = A(i,:);
    if i==1
        plot(x,y,'-|',...
            'LineWidth',2,...
            'Color', color1,...
            'MarkerSize',2,... 
            'MarkerFaceColor',color1);
    elseif i==2
        plot(x,y,'-|',...
            'LineWidth',2,...
            'Color', color2,...
            'MarkerSize',2,... 
            'MarkerFaceColor',color2);
    else
        plot(x,y,'-|',...
            'LineWidth',2,...
            'Color', color3,...
            'MarkerSize',2,... 
            'MarkerFaceColor',color3);
    end
    hold on;

end
title(['Bi=' num2str(pup)], 'fontsize', 12);
xlabel('x, м', 'fontsize', 12);
ylabel('T, K', 'fontsize', 12);
xlim([-0.01*x(101),x(101)+0.01*x(101)]);
ylim([400,800]);
legend('t=0.1∙Fo_m_a_x', 't=0.5∙Fo_m_a_x', 't=0.9∙Fo_m_a_x','Location','northeast');
hold off;
saveas(v, [folders 'temp_x_alpha' int2str(lop) '.png']);
close
end



FP = load(['res/fo_points.txt']);
fnum = FP(:,1);
fpoints = FP(:,2);
fzero = FP(:,3);
for j=1:2:3
k=figure;
lop=round(AlphaCoord(j));
pup=BiCoord(j);

x=rand(101,1);
y=rand(101,1);

for i=1:1:3

    C = load([folder1 'alpha' int2str(lop) '_temp_t.txt']);
    D = load([folder1 'alpha' int2str(lop) '_time.txt']);

    x(1) = D(1,1);
    y(1) = C(1,i);
    x(101) = D(time_lim(j),1);
    y(101) = C(time_lim(j),i);
    ii_t=1;
    for ii=2:100
        ii_t=round(ii_t+time_lim(j)/99);
        x(ii) = D(ii_t,1);
        y(ii) = C(ii_t,i);
    end

    if i==1
            plot(x,y,'-|',...
                'LineWidth',2,...
                'Color', color1,...
                'MarkerSize',2,... 
                'MarkerFaceColor',color1);
    elseif i==2
            plot(x,y,'-|',...
                'LineWidth',2,...
                'Color', color2,...
                'MarkerSize',2,... 
                'MarkerFaceColor',color2);
    else
            plot(x,y,'-|',...
                'LineWidth',2,...
                'Color', color3,...
                'MarkerSize',2,... 
                'MarkerFaceColor',color3);
    end
    hold on;
end

plot(fpoints,fzero,'|', ... 
    'MarkerSize',10);
hold on;
title(['Bi=' num2str(pup)], 'fontsize', 12);
xlabel('t, c', 'fontsize', 12);
ylabel('T, K', 'fontsize', 12);
xlim([-0.01*x(101),x(101)+0.01*x(101)]);
ylim([400,800]);
legend('x=0', 'x=d/4', 'x=d/2','Location','northeast');
hold off;
saveas(k, [folders 'temp_t_alpha' int2str(lop) '.png']);
close
end