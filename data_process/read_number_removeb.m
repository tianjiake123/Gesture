%����Ƶԭʼ�ź�ͨ��FFTת����Ƶ����ȡ���롢�ٶ���Ϣ��������ͨ����ȡ�Ƕ���Ϣ�����Ŀ���ķ�λ��Ϣ
%��֡��λ��Ϣ��������p_file�ļ���
%ר�ø������ʣ�
%       ɢ��㣺һ��Ŀ��һ��ᱻ�״���ɼ���ɢ�����ģ����ǳ�֮Ϊɢ���
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

p_file= dir('.\raw_data\number\*.dat');    %������������ԭʼ����(raw data)�ļ���������������dat�ļ�����C++�ļ���ȡ��õģ�bin�ļ�����mmwave studio����ɼ���õģ����߱���������
file_name = {p_file.name};           %��¼����·��������
file_folder = {p_file.folder};

load('b.mat', 'background')           %����Ԥ��ı���RDM�����ڱ�������������read_background.m�����ǽ��ɼ����ı���(Ŀ�겻��̽�ⷶΧ������ת��ΪRDM���洢

for ko = 1:length(file_name)          %����raw data�ļ��У����������ļ�
    fid = fopen([file_folder{ko} '\' file_name{ko}],'r');%adc_microdoppler
    det_matrix = zeros(64,128);
    for ki = 1:10000                   %һ���ļ�������֡���ݣ���֡����
        data = fread(fid,[2*Num_rx*Num_samples,Num_chirps],'int16');    %��ȡԭʼadc���ݣ�ÿ�ζ�ȡһ֡
        [a, b] = size(data);           %һ֡���ݲ�����������
        if isempty(data) || ~(a==1024) || ~(b==128)
            break
        end
        %һ֡adc���ݵķֲ���chirp1-chirp2-chirp3-...��ÿ��chirp�а����������ߺͲ���������ݣ�
        %��һ��chirp�����ݴ�СΪN_sample��N_rx��N_tx��2(I/Qֵ��
        %���е������Ų�Ϊrx0i0 rx0i1 rx0q0 rx0q1 rx0i2 rx0i3 rx0q2 rx0q3...
        %�ĸ�һ����и�ֵ
        %����������Ϊһ����ά�������飺range*angle*velocity
        data_c(1:2:end,:) = data(1:4:end,:) + 1j * data(3:4:end,:);
        data_c(2:2:end,:) = data(2:4:end,:) + 1j * data(4:4:end,:);
        data_adc = reshape(data_c,[Num_samples,Num_rx,Num_chirps]);
        %��ÿ������ͨ����2D-FFT�����ÿ��ͨ����RDM
        for ks = 1:Num_rx
            aa = squeeze(data_adc(:,ks,1:end));
            bb = fft2(aa.*range_window.*velocity_window,fft_range,fft_velocity);
            cc(:,:,ks) = fftshift(bb,2);                                            %�����ٶ���������fftshiftʵ����һ�����ݵĽ���λ�ã����ٶ�Ϊ0�Ƶ����м�
        end  
        %ȡһ��ͨ�������ݣ�����FFT�����Ŀ�������ֵ������һ����˵����ͬͨ���ķ�ֵ���࣬����һ��ֻ��һ��ͨ����RDM��������ֵ
        aa = squeeze(data_adc(:,1,1:end));
        data_fft = fft2(aa.*range_window.*velocity_window,fft_range,fft_velocity);
        da_fft = fftshift(data_fft,2);da_fft_mod = abs(da_fft);
        %��������������ǰRDM��ȥԤ��ı���RDM
        rdm = da_fft_mod - background;    
        %������ֹ����ĵ��ֵ������0
        rdm(Range_set:end,:) = 0;
        
        %��ֵ����ֵ�������ҵ������㣬һ����˵��RDM�Ϸ�ֵ��ֵ���Ӧ��ɢ��㣬���ͨ����ֵ�������ҵ�ɢ���
        det_matrix(1:Range_set,:)= isjizhi(rdm(1:Range_set,:)).*rdm(1:Range_set,:);
        det_matrix = det_matrix.*(det_matrix>Detect_noise);                     %��ֵ���������ֵ�������������ţ��˴���ֵ����Ϊ����ֵ
        det_matrix(:,65) = det_matrix(:,65).*(det_matrix(:,65)>Static_noise);   %65=128/2,���65��Ӧ�ٶ�Ϊ0������ֹ���ĵ㣬������Ӧһ���ϴ�ķ�ֵ�����������һ���ϴ����ֵ���ų�����ţ�Ч����̫���ԣ����Ż�
        det_matrix(1:6,65) = 0;                                                 %��ֹ���������ֵ��0����ֹ���Ÿ���
        det_matrix(11:end,65) = 0;                                              %�Ծ�ֹĿ�����һ����С�Ľ�ֹ���룬��ֹ����ĸ���
        
        %ʵ���ϣ�����Ӱ��̽��������64��66���е����ݣ��������˶����壬���ֱ�ӽ���̬���弴65����0��
        %64��66��Ϊ�߽�㣬�Ӷ��޷���Ϊ��ֵ�㣬����Ŀ���������˶�ʱ�������ٶ�̫���������޷�����⵽����5��7����һ��
        %���ֻ��63-67����������������ֵ���������һ��

        [ind_x,ind_y,~] =find(det_matrix > Detect_noise);                       %��ȡ��Ӧ����ֵ��x��Ӧ���룬y��Ӧ�ٶ�
        ind_v = ind_y - 65;                                                     %�����ٶ�
        
        if length(ind_y)>1  % Ҫ��֡�Ĳ���Ŀ��ɢ�����������������ǿ�ȶ���
            %���ֵ���ڷ�ֵ��Ȩ�����ģ���Ϊ��֡��̽����
            temp_v = 0;temp_r = 0;temp_mod = 0;temp_az = 0;temp_el = 0;temp_el1 = 0;temp_ver_vel = 0;
            for as = 1:length(ind_y)
                temp_v = det_matrix(ind_x(as),ind_y(as)) * ind_v(as) + temp_v;
                temp_r = det_matrix(ind_x(as),ind_y(as)) * ind_x(as) + temp_r;
                temp_mod = det_matrix(ind_x(as),ind_y(as)) + temp_mod;
                jiaodu_az = angle(cc(ind_x(as),ind_y(as),4)) - angle(cc(ind_x(as),ind_y(as),1));     %����λ��Ƕȣ���������ͨ���Ĳ����й�
                                                                                                     %����Ŀʹ�õ�IWR6842 ods�״�
                                                                                                     % MIMO�����Ų�Ϊ   R1 R4    T1 T2
                                                                                                     %                 R2 R3     
                                                                                                     %��ЧΪ           1 4 5 8
                                                                                                     %                 2 3 6 7
                                                                                                     % ��˿�ȡ1 4 5 8ͨ������AZ  
                jiaodu_el = angle(cc(ind_x(as),ind_y(as),2)) - angle(cc(ind_x(as),ind_y(as),1)) - pi;%- pi%isk��Ҫ����ͨ��������ĿĿǰ��û�õ���λ��
               
                %�޶��Ƕ���-pi��pi֮��
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
            if temp_mod > Sum_noise             %�����յķ�ֵ�ʹ�����ֵʱ����֡��Ч����¼Ŀ����Ϣ
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










