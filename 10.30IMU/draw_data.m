function draw_data(X1, Y1, X2, Y2)
%CREATEFIGURE(X1, Y1, X2, Y2)
%  X1:  plot x 数据的向量
%  Y1:  plot y 数据的向量
%  X2:  plot x 数据的向量
%  Y2:  plot y 数据的向量

%  由 MATLAB 于 15-Nov-2023 20:03:05 自动生成

% 创建 figure
figure1 = figure('WindowState','maximized');

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 创建 plot
plot(X1,Y1,'DisplayName','纯惯导解算结果','LineWidth',1,...
    'Color',[0.850980392156863 0.325490196078431 0.0980392156862745]);

% 创建 plot
plot(X2,Y2,'DisplayName','参考真值','LineStyle','--',...
    'Color',[0.301960784313725 0.745098039215686 0.933333333333333]);

% 创建 ylabel
ylabel({'N(m)'});

% 创建 xlabel
xlabel({'E(m)'});

% 创建 title
title({'实验数据纯惯导解算结果和参考结果'});

% 取消以下行的注释以保留坐标区的 X 范围
xlim(axes1,[-25.6645379162378 38.7759942903615]);
% 取消以下行的注释以保留坐标区的 Y 范围
ylim(axes1,[-14.5301341023262 23.0601763515234]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% 创建 legend
legend(axes1,'show');
end
