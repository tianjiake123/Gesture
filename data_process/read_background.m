%计算背景RDM，并保存
clc;clear;close all;
% 雷达波形参数
Num_rx = 8;             %接收天线数，与硬件配置一致，本项目的实际接收天线数为4，但是因为采数时设置了2T4R MIMO技术，等同为8个接收天线，但本项目中算法实际只使用了2根接收天线
Num_samples = 64;       %采样点数，与目标径向距离有关
Num_chirps = 128;       %chirp数，此参数与目标径向速度分辨率有关
fft_range = 64;         %距离维FFT点数，一般与采样数一致，不一致时，相当于进行了插值或下采样，一般是2的幂次。
fft_velocity = 128;     %速度维FFT点数，一般与chirp数一致，不一致时，相当于进行了插值或下采样，一般是2的幂次。

Range_set = 20;         %动作截止距离，目前距离分辨率为5.4cm，因此最大距离为5.4×20=108cm，这里不认为手的移动距离会超过这个范围
Detect_noise = 1e4;     %判断目标点是否有效的幅值阈值（散射动态点）
Static_noise = 1e5;     %判断目标点是否有效的幅值阈值（散射静态点）
Sum_noise = 1e5;        %判断目标点是否有效的幅值阈值（质心点）
Num_flag = 5;           %记录的目标点信息数量，此处每帧目标点保存信息包括：（径向）距离、（径向）速度、角度、垂直距离、垂直速度

data_c = zeros(Num_rx*Num_samples,Num_chirps);                      %装填一帧cube数据用
range_window = hamming(Num_samples)*ones(1,fft_velocity);           %加窗，此处选用hamming窗，后续FFT用，最好不要用矩形窗（即直接截取）
velocity_window = ones(Num_samples,1)*hamming(fft_velocity)';

fid = fopen('.\raw_data\b_ods_high.dat','r');                 %读取背景原始文件
det_matrix = zeros(64,128);
for ki = 1:10000                                                    %逐帧处理，文件读完会跳出循环，所以直接设置了一个较大值
    data = fread(fid,[2*Num_rx*Num_samples,Num_chirps],'int16');    %读取一帧原始数据
    [a, b] = size(data);                                            %一帧数据不完整，舍弃
    if isempty(data) || ~(a==1024) || ~(b==128)
        break
    end
    data_c(1:2:end,:) = data(1:4:end,:) + 1j * data(3:4:end,:);
    data_c(2:2:end,:) = data(2:4:end,:) + 1j * data(4:4:end,:);
    data_adc = reshape(data_c,[Num_samples,Num_rx,Num_chirps]);%range*angle*velocity
    for ks = 1
        aa = squeeze(data_adc(:,ks,1:end));
        data_fft = fft2(aa.*(hamming(Num_samples)*ones(1,fft_velocity)).*velocity_window,fft_range,fft_velocity);
        da_fft = fftshift(data_fft,2);da_fft_mod = abs(da_fft);
        da_fft_mod(Range_set:end,:) = 0;                        %超出范围的区域置0
%         计算观察RDM
%         mesh(da_fft_mod(1:Range_set,:))
        rdm(:,:,ks)= da_fft_mod;
    end 
    background(:,:,ki) = mean(rdm,3);                            %多通道求均值    
end
background = mean(background, 3);                                %多帧求均值
%雷达信号有不稳定性，积累均值可以使获得的RDM更稳定，接近真实值
save('b.mat', 'background');









