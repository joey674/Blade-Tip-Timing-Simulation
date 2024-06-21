clear;

% 加载数据
load('DL/peak.mat', 'P', 'T');

% 打印加载的数据维度
disp('加载数据维度:');
disp(['P: ', num2str(size(P))]);
disp(['T: ', num2str(size(T))]);

% 训练线性神经网络
lr = maxlinlr(P); % 获取最大学习速率
net = newlin(minmax(P), 1, 0, lr); % 建立线性神经网络
net.trainParam.epochs = 1000; % 训练：做多500次
net.trainParam.goal = 0.04; % 训练误差设定为0.04
net = train(net, P, T); % 训练
Y = sim(net, P); % 仿真

% 显示仿真结果
disp('仿真结果:');
% disp(Y);

% 计算并显示训练误差
train_error = perform(net, T, Y);
disp(['训练误差: ', num2str(train_error)]);

% 绘制实际标签和仿真结果
num_samples = size(P, 1);

for i = 1:num_samples
    figure('Name', ['Sample ' num2str(i)], 'NumberTitle', 'off');
    subplot(3, 1, 1);
    plot(P(i, :), 'b');
    title(' Input');
    xlabel('Index');
    ylabel('Magnitude');
    
    subplot(3, 1, 2);
    plot(T(i, :), 'ro', 'DisplayName', 'Actual Peaks');

    subplot(3, 1, 3);
    plot(Y(i, :), 'go', 'DisplayName', 'Predicted Peaks');

    legend;
    hold off;
end
