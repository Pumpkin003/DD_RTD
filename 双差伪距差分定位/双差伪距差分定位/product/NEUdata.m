clear;
clc;

% 导入数据
data = importdata("3RecCoorDataSTR2.txt");

% 计算数据的平均值和标准差
meandata = mean(data);
maxdata = max(data);

% 历元时间间隔（秒）
epoch_interval = 30;

% 生成时间轴（以秒为单位）
time_seconds = (0:size(data,1)-1) * epoch_interval;

% 将时间转换为日期时间格式
time = datetime(time_seconds, 'ConvertFrom', 'epochtime', 'Epoch', '2024-03-23', 'Format', 'HH:mm:ss');

% 创建一个新的图形窗口
figure('Color', 'w'); % 设置背景颜色为白色
hold on

% 绘制在N、E、U方向的误差
plot(time, data(:,1), '-b', 'LineWidth', 0.1, 'DisplayName', '在N方向误差');
plot(time, data(:,2), '-g', 'LineWidth', 0.1, 'DisplayName', '在E方向误差');
plot(time, data(:,3), '-r', 'LineWidth', 0.1, 'DisplayName', '在U方向误差');

% 添加均值水平线
x_limits = [time(1) time(end)]; % 获取当前x轴的范围
plot(x_limits, [meandata(1) meandata(1)], '--b', 'LineWidth', 2, 'DisplayName', sprintf('在N方向误差均值 (%.2f)', meandata(1)));
plot(x_limits, [meandata(2) meandata(2)], '--g', 'LineWidth', 2, 'DisplayName', sprintf('在E方向误差均值 (%.2f)', meandata(2)));
plot(x_limits, [meandata(3) meandata(3)], '--r', 'LineWidth', 2, 'DisplayName', sprintf('在U方向误差均值 (%.2f)', meandata(3)));

% 添加标题和轴标签
title('定位误差分析', 'FontSize', 14, 'FontWeight', 'bold')
xlabel('时间 (HH:mm:ss)', 'FontSize', 12)
ylabel('定位误差 (米)', 'FontSize', 12)

% 添加网格
grid on
grid minor

% 添加图例
legend('show', 'Location', 'best')

% 设置轴的字体大小
set(gca, 'FontSize', 10)

% 关闭保持状态
hold off
