%将轨迹图像转换为像素图片
load('p_file.mat', 'p_file');
%像素图片的尺寸，这里选择28是为了和手写数字数据集MNIST保持一致
w = 28;
h = 28;

% 边界值，此处是经验值
x_min = -5;
x_max = 5;
y_min = 4;
y_max = 14;
%
for j = 1:length(p_file)

    data = p_file(j).data;
    % 读取存储的目标轨迹
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
    img = zeros(28,28);
    %轨迹方位值转化未像素点索引值
    x = round(((s_hor_range-x_min)/(x_max-x_min))*28);
    y = round((1-(s_ver_range-y_min)/(y_max-y_min))*28);

    %对应像素点置1
    for i = 1:length(x)
        img(x(i),y(i)) = 1;
    end
    %记录转换后的值,存成dat比较容易给c++程序读取，用imwrite(img,'result.jpg');可以存成图片
    if ~exist('./imgs')
        mkdir('./imgs')
    end
    fid=fopen(['./imgs/' num2str(p_file(j).label) '.dat'],'wb');
    fwrite(fid,single(img),'float');
    fclose(fid);
    imshow(img')
end

