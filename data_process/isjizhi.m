function [outputArg1] = isjizhi(inputArg1)
%isjizhi ��ֵ����
%   inputArg1��һ����ά���飬��������������������ֵ��ֵ������������ֵ��ע��߽�㲻�ᱻ��Ϊ��ֵ��
    aa1 = islocalmax(inputArg1,1);
    aa2 = islocalmax(inputArg1,2);
    outputArg1 = aa1.*aa2;
end

