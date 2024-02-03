%将中频原始信号通过FFT转换到频域，提取距离、速度信息，多天线通道提取角度信息，输出目标点的方位信息
%多帧方位信息被储存在p_file文件中
%专用概念名词：
%       散射点：一个目标一般会被雷达检测成几个散射中心，我们称之为散射点
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

p_file= dir('.\raw_data\number\*.dat');    %查找所有离线原始数据(raw data)文件，用于批量处理，dat文件是用C++文件读取获得的，bin文件是用mmwave studio程序采集获得的，两者本质无区别
file_name = {p_file.name};           %记录数据路径和名称
file_folder = {p_file.folder};

load('b.mat', 'background')           %加载预存的背景RDM，用于背景减除，运行read_background.m，我们将采集到的背景(目标不在探测范围）数据转换为RDM并存储

for ko = 1:length(file_name)          %遍历raw data文件夹，挨个处理文件
    fid = fopen([file_folder{ko} '\' file_name{ko}],'r');%adc_microdoppler
    det_matrix = zeros(64,128);
    for ki = 1:10000                   %一个文件包括多帧数据，逐帧处理
        data = fread(fid,[2*Num_rx*Num_samples,Num_chirps],'int16');    %读取原始adc数据，每次读取一帧
        [a, b] = size(data);           %一帧数据不完整，舍弃
        if isempty(data) || ~(a==1024) || ~(b==128)
            break
        end
        %一帧adc数据的分布：chirp1-chirp2-chirp3-...，每个chirp中包括所有天线和采样点的数据，
        %即一个chirp的数据大小为N_sample×N_rx×N_tx×2(I/Q值）
        %其中的数据排布为rx0i0 rx0i1 rx0q0 rx0q1 rx0i2 rx0i3 rx0q2 rx0q3...
        %四个一组进行赋值
        %将数据重排为一个三维复数数组：range*angle*velocity
        data_c(1:2:end,:) = data(1:4:end,:) + 1j * data(3:4:end,:);
        data_c(2:2:end,:) = data(2:4:end,:) + 1j * data(4:4:end,:);
        data_adc = reshape(data_c,[Num_samples,Num_rx,Num_chirps]);
        %对每个天线通道做2D-FFT，获得每个通道的RDM
        for ks = 1:Num_rx
            aa = squeeze(data_adc(:,ks,1:end));
            bb = fft2(aa.*range_window.*velocity_window,fft_range,fft_velocity);
            cc(:,:,ks) = fftshift(bb,2);                                            %由于速度有正负，fftshift实际是一个数据的交换位置，将速度为0移到了中间
        end  
        %取一个通道的数据，计算FFT，获得目标的索引值，由于一般来说，不同通道的幅值相差不多，我们一般只对一个通道的RDM查找索引值
        aa = squeeze(data_adc(:,1,1:end));
        data_fft = fft2(aa.*range_window.*velocity_window,fft_range,fft_velocity);
        da_fft = fftshift(data_fft,2);da_fft_mod = abs(da_fft);
        %背景减除，将当前RDM减去预存的背景RDM
        rdm = da_fft_mod - background;    
        %超过截止距离的点幅值将被置0
        rdm(Range_set:end,:) = 0;
        
        %极值（峰值）搜索找到动检测点，一般来说，RDM上幅值极值点对应了散射点，因此通过峰值搜索查找到散射点
        det_matrix(1:Range_set,:)= isjizhi(rdm(1:Range_set,:)).*rdm(1:Range_set,:);
        det_matrix = det_matrix.*(det_matrix>Detect_noise);                     %幅值必须大于阈值，减少噪声干扰，此处阈值设置为经验值
        det_matrix(:,65) = det_matrix(:,65).*(det_matrix(:,65)>Static_noise);   %65=128/2,因此65对应速度为0（即静止）的点，人体会对应一个较大的幅值，因此设置了一个较大的阈值以排除其干扰，效果不太明显，待优化
        det_matrix(1:6,65) = 0;                                                 %静止物体近处幅值置0，防止串扰干扰
        det_matrix(11:end,65) = 0;                                              %对静止目标会有一个更小的截止距离，防止人体的干扰
        
        %实际上，真正影响探测结果的是64、66两列的数据，即低速运动物体，如果直接将静态物体即65处置0，
        %64、66成为边界点，从而无法成为极值点，导致目标做切向运动时，径向速度太低以至于无法被检测到（如5、7）的一横
        %因此只对63-67处做背景减除并极值搜索，结果一致

        [ind_x,ind_y,~] =find(det_matrix > Detect_noise);                       %获取对应索引值，x对应距离，y对应速度
        ind_v = ind_y - 65;                                                     %正负速度
        
        if length(ind_y)>1  % 要求单帧的测量目标散射点至少有两个，增强稳定性
            %获得值关于幅值加权求质心，作为此帧的探测结果
            temp_v = 0;temp_r = 0;temp_mod = 0;temp_az = 0;temp_el = 0;temp_el1 = 0;temp_ver_vel = 0;
            for as = 1:length(ind_y)
                temp_v = det_matrix(ind_x(as),ind_y(as)) * ind_v(as) + temp_v;
                temp_r = det_matrix(ind_x(as),ind_y(as)) * ind_x(as) + temp_r;
                temp_mod = det_matrix(ind_x(as),ind_y(as)) + temp_mod;
                jiaodu_az = angle(cc(ind_x(as),ind_y(as),4)) - angle(cc(ind_x(as),ind_y(as),1));     %求相位差（角度），与天线通道的布局有关
                                                                                                     %本项目使用的IWR6842 ods雷达
                                                                                                     % MIMO天线排布为   R1 R4    T1 T2
                                                                                                     %                 R2 R3     
                                                                                                     %等效为           1 4 5 8
                                                                                                     %                 2 3 6 7
                                                                                                     % 因此可取1 4 5 8通道计算AZ  
                jiaodu_el = angle(cc(ind_x(as),ind_y(as),2)) - angle(cc(ind_x(as),ind_y(as),1)) - pi;%- pi%isk需要调整通道，本项目目前还没用到方位角
               
                %限定角度在-pi到pi之间
                if jiaodu_az > pi
                    jiaodu_az = jiaodu_az - 2*pi;
                end

                if jiaodu_az < -pi
                    jiaodu_az = jiaodu_az + 2*pi;
                end

                if jiaodu_el > pi
                    jiaodu_el = jiaodu_el - 2*pi;
                end

                if jiaodu_el < -pi
                    jiaodu_el = jiaodu_el + 2*pi;
                end
                
                temp_az = det_matrix(ind_x(as),ind_y(as)) * jiaodu_az + temp_az;
                temp_el = det_matrix(ind_x(as),ind_y(as)) * jiaodu_el + temp_el;
                temp_ver_vel = det_matrix(ind_x(as),ind_y(as))*ind_v(as).*sqrt(1-(jiaodu_az/pi).^2) + temp_ver_vel;  
            end
            temp_v = temp_v/temp_mod;
            temp_r = temp_r/temp_mod;
            temp_az = temp_az/temp_mod;
            temp_el = temp_el/temp_mod;
            temp_el1 = temp_el1/temp_mod;
            ver_range = temp_r.*sqrt(1-(temp_az/pi).^2); 
            ver_vel = temp_ver_vel/temp_mod; 
            if temp_mod > Sum_noise             %当最终的幅值和大于阈值时，此帧有效，记录目标信息
                flag(ki,1:Num_flag) = [temp_v,temp_r,temp_az, ver_range, ver_vel];
            else
                flag(ki,1:Num_flag) = 0;
            end
        else
            flag(ki,1:Num_flag) = 0;
        end
    end
    fclose(fid);
    p_file(ko).data =  flag;
    p_file(ko).label = str2num(p_file(ko).name(1));

end
save('p_file.mat', 'p_file');










