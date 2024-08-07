% 读文件 & 变量赋值
Flag = 1;                               % flag = 1表示灰色。flag = 2表示灰绿
ArrayLength = 8;                        % PossibleRootNodes/PossibleRootpipes数组长度
FirstRootNodes = 41;                    % PossibleRootNodes起始元素
FirstRootpipes = 58;                    % PossibleRootpipes起始元素
NumberOfOutlets = 2;                    % 控制出口数量，修改即可代表有几个出口, 一定要小于等于ArrayLength

Outletindex = zeros(1, ArrayLength); 
% 随机选择x个索引  
indices = unique(randi(ArrayLength, NumberOfOutlets, 1));  
% 将选定的索引位置上的元素赋值为1  
Outletindex(indices) = 1;  

%LIDchushizhongqun=randi([0,1],400,202);

load("chushihzhongqun.mat")
Population=LIDchushizhongqun;
GenomeLength = 110;
%GenomeLength = 202;
MaxGenerations = 5;
PopulationSize = 400;
options = optimoptions('gamultiobj','MutationFcn',@mutationadaptfeasible,'MaxGenerations',MaxGenerations,'InitialPopulationMatrix',Population,...
    'PopulationSize',PopulationSize,'PlotFcn',{@gaplotscorediversity,@gaplotspread,@gaplotpareto},'PopulationType','bitstring');

% BIM=dlmread('BIM.txt');                 % 管道矩阵
% loops=dlmread('Loops.txt');             % 集水区矩阵
% condults=dlmread('[CONDUITS].txt');     % 管道上下游节点矩阵
% nodescordinates=dlmread('nodescordinates.txt');
% Commercial=dlmread('Commercial.txt');
% InputfileOrg = regexp( fileread('AhvazNull.inp'), '\n', 'split');
% Initial=dlmread('[SUBCATCHMENTS].txt')
% ToSearch=regexp( fileread('[SUBCATCHMENTS] - T.txt'), '\n', 'split');
% SUBAREAS=dlmread('[SUBAREAS].txt');
% ToSearch=regexp( fileread('[SUBAREAS]-T.txt'), '\n', 'split');
% INFILTRATION=dlmread('[INFILTRATION].txt');
% ToSearch=regexp( fileread('[INFILTRATION]-T.txt'), '\n', 'split');
% Polygons=dlmread('[Polygons].txt');
% ToSearch=regexp( fileread('[Polygons]-T.txt'), '\n', 'split');
% ToSearch=regexp( fileread('[JUNCTIONS]-T.txt'), '\n', 'split');
% ToSearch=regexp( fileread('[OUTFALLS]-T.txt'), '\n', 'split');
% ToSearch=regexp( fileread('[CONDUITS]-T.txt'), '\n', 'split');
% ToSearch=regexp( fileread('[XSECTIONS]-T.txt'), '\n', 'split');
% fid = fopen('LIMatF2.inp', 'w');  %这个.inp文件是什么？
% NodeEle=dlmread('NodeEle.txt');

fitnessFunction = @(x)CANew_hy(x, Outletindex, ArrayLength, FirstRootNodes, FirstRootpipes, Flag);  

if Flag == 1
    [inp,fval,exitflag,output,Population,scores]=gamultiobj(fitnessFunction,GenomeLength,[],[],[],[],[],[],options); %入参
    % 存储最优解的变量值
    save('inp.mat', 'inp');
else
    [dinp,fval,exitflag,output,Population,scores]=gamultiobj(fitnessFunction,GenomeLength,[],[],[],[],[],[],options);
    save('dinp.mat', 'dinp');
end


