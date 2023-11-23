%% 从文本文件中导入数据
% 用于从以下文本文件中导入数据的脚本:
%
%    filename: E:\文件\复习\惯性导航\IMU10.8\duizhun.ASC
%
% 由 MATLAB 于 2023-10-13 15:27:34 自动生成

%% 设置导入选项并导入数据
function[data]=readfile(filepath)
    opts = delimitedTextImportOptions("NumVariables", 12);
    
    % 指定范围和分隔符
    opts.DataLines = [1, Inf];
    opts.Delimiter = ["*", ","];
    
    % 指定列名称和类型
    opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "Var12"];
    opts.SelectedVariableNames = ["VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11"];
    opts.VariableTypes = ["string", "string", "string", "string", "string", "double", "double", "double", "double", "double", "double", "string"];
    
    % 指定文件级属性
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % 指定变量属性
    opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var12"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var12"], "EmptyFieldRule", "auto");

    % 导入数据
    data_value = readtable(filepath, opts);
    %% 转换为输出类型
    data=table2array(data_value);
    %% 清除临时变量
    clear opts
end
