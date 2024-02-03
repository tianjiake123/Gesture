function output = IIR(input, factor)
%IIR 递归函数，平滑滤波
%   input：需要平滑的序列
%   factor：保留系数
output = zeros(length(input),1);
temp = input(1);
output(1) = temp;
for i = 2:length(input)
    temp = temp * factor + input(i) *(1-factor);
    output(i) = temp;
end
end

