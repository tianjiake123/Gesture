%绘制保存的p_file文件中的目标轨迹图像，检查结果
clc;clear;close all;
load('p_file.mat', 'p_file');

for j = 1:length(p_file)

    data = p_file(j).data;
    %读取存储的目标轨迹
    vel = data(:,1);
    range = data(:,2);
    phase = data(:,3);
    ver_range = data(:,4);
    ver_vel = data(:,5);
    
    %用vel为标准要求有探测到动态物体
    phase(vel==0) = [];
    ver_range(vel==0)=[];
    ver_vel(vel==0)=[];
    range(vel==0)=[];
    vel(vel==0) = [];
    
    %计算水平位置
    hor_range = -range.*phase/pi;
    %平滑轨迹
    s_hor_range = IIR(hor_range,0.8);
    s_ver_range = IIR(ver_range,0.8);
    %绘制图像
    subplot(2,5,j)
    %平滑曲线
    plot(s_hor_range, s_ver_range, 'k', 'linewidth', 2)
    hold on
    %未平滑的点
    plot(hor_range, ver_range,'bo','MarkerFaceColor', 'b', 'MarkerSize', 3)

    title( num2str(p_file(j).label))

end