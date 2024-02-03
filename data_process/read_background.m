%���㱳��RDM��������
clc;clear;close all;
% �״ﲨ�β���
Num_rx = 8;             %��������������Ӳ������һ�£�����Ŀ��ʵ�ʽ���������Ϊ4��������Ϊ����ʱ������2T4R MIMO��������ͬΪ8���������ߣ�������Ŀ���㷨ʵ��ֻʹ����2����������
Num_samples = 64;       %������������Ŀ�꾶������й�
Num_chirps = 128;       %chirp�����˲�����Ŀ�꾶���ٶȷֱ����й�
fft_range = 64;         %����άFFT������һ���������һ�£���һ��ʱ���൱�ڽ����˲�ֵ���²�����һ����2���ݴΡ�
fft_velocity = 128;     %�ٶ�άFFT������һ����chirp��һ�£���һ��ʱ���൱�ڽ����˲�ֵ���²�����һ����2���ݴΡ�

Range_set = 20;         %������ֹ���룬Ŀǰ����ֱ���Ϊ5.4cm�����������Ϊ5.4��20=108cm�����ﲻ��Ϊ�ֵ��ƶ�����ᳬ�������Χ
Detect_noise = 1e4;     %�ж�Ŀ����Ƿ���Ч�ķ�ֵ��ֵ��ɢ�䶯̬�㣩
Static_noise = 1e5;     %�ж�Ŀ����Ƿ���Ч�ķ�ֵ��ֵ��ɢ�侲̬�㣩
Sum_noise = 1e5;        %�ж�Ŀ����Ƿ���Ч�ķ�ֵ��ֵ�����ĵ㣩
Num_flag = 5;           %��¼��Ŀ�����Ϣ�������˴�ÿ֡Ŀ��㱣����Ϣ�����������򣩾��롢�������ٶȡ��Ƕȡ���ֱ���롢��ֱ�ٶ�

data_c = zeros(Num_rx*Num_samples,Num_chirps);                      %װ��һ֡cube������
range_window = hamming(Num_samples)*ones(1,fft_velocity);           %�Ӵ����˴�ѡ��hamming��������FFT�ã���ò�Ҫ�þ��δ�����ֱ�ӽ�ȡ��
velocity_window = ones(Num_samples,1)*hamming(fft_velocity)';

fid = fopen('.\raw_data\b_ods_high.dat','r');                 %��ȡ����ԭʼ�ļ�
det_matrix = zeros(64,128);
for ki = 1:10000                                                    %��֡�����ļ����������ѭ��������ֱ��������һ���ϴ�ֵ
    data = fread(fid,[2*Num_rx*Num_samples,Num_chirps],'int16');    %��ȡһ֡ԭʼ����
    [a, b] = size(data);                                            %һ֡���ݲ�����������
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
        da_fft_mod(Range_set:end,:) = 0;                        %������Χ��������0
%         ����۲�RDM
%         mesh(da_fft_mod(1:Range_set,:))
        rdm(:,:,ks)= da_fft_mod;
    end 
    background(:,:,ki) = mean(rdm,3);                            %��ͨ�����ֵ    
end
background = mean(background, 3);                                %��֡���ֵ
%�״��ź��в��ȶ��ԣ����۾�ֵ����ʹ��õ�RDM���ȶ����ӽ���ʵֵ
save('b.mat', 'background');









