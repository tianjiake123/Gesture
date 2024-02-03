function [outputArg1] = isjizhi(inputArg1)
%isjizhi 极值搜索
%   inputArg1是一个二维数组，在行列两个方向搜索幅值极值，并返回索引值，注意边界点不会被视为极值点
    aa1 = islocalmax(inputArg1,1);
    aa2 = islocalmax(inputArg1,2);
    outputArg1 = aa1.*aa2;
end

