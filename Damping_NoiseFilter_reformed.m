

function magn_modi = Damping_NoiseFilter_reformed(magn,peaks_idx,weights_idx)   
    magn = smoothdata(magn,'movmean',2); % 用一点点的mov去消除过于难以接受的噪音;
    magn_gauss = smoothdata(magn, 'gaussian', 100);% 我们会使用的模型,但是尺度有问题,底部保留,峰值处还可以再高;
    
    adjust_conf = 4/6;% 这里设置成0就是对整个范围拉伸;这个值代表拉伸范围
    max_adjustRatio = 0.5;
    % 初始化修改后的数据为高斯平滑数据
    magn_modi = magn_gauss;  
    % 遍历每个找到的峰值
    for j = 1:length(peaks_idx)
        peakCenter = peaks_idx(j);
        startIndex = round( weights_idx(j,1) + (peaks_idx(j)-weights_idx(j,1))*adjust_conf );
        endIndex   = round( weights_idx(j,2) - (weights_idx(j,2)-peaks_idx(j))*adjust_conf );
        % 获取峰值区间内的数据
        peakData = magn_gauss(startIndex:endIndex);
        % 计算峰值区间内的最大值
        peakMax = max(peakData);
        % 对峰值区间内的数据进行线性减弱的尺度调整
        for i = startIndex:endIndex
            % 计算当前点到峰值中心的距离比例
            distanceRatio = abs(i - peakCenter) / (endIndex - startIndex);
            % 计算调整比例，越靠近峰值中心，调整比例越大
            adjustRatio = 1 - distanceRatio; % 距离中心越近，adjustRatio越大
            if adjustRatio > max_adjustRatio
                adjustRatio = max_adjustRatio;
            end
            % 调整当前点的尺度
            magn_modi(i) = magn_modi(i) / peakMax * max(magn(startIndex:endIndex)) * adjustRatio + magn_modi(i) * (1 - adjustRatio);
        end
    end
end
