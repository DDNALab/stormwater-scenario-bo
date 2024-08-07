function obj=CANew_hy(inp, outletindex, arrayLength, firstRootNodes, firstRootpipes, flag)
  %%%%%%%
dinp = inp;  % 仅在灰绿使用
    if flag ~= 1
       try
            inp = load("inp.mat");
            inp = inp.inp;
       catch
            disp('灰色inp文件缺失，请添加灰绿inp文件后重试');  % 处理加载失败的情况
       end
    end


    % if numberOfOutlets==1 
    % outletindex=[0 0 0 0 0 1 0 0]; 
    % elseif numberOfOutlets==2
    % outletindex=[0 0 1 0 0 1 0 0];
    % elseif numberOfOutlets==3
    % outletindex=[1 0 1 0 0 1 0 0];
    % elseif numberOfOutlets==4        
    % outletindex=[1 0 1 0 1 0 1 0];
    % elseif numberOfOutlets==5
    % outletindex=[1 1 0 1 0 1 0 1];
    % elseif numberOfOutlets==6
    % outletindex=[1 1 0 1 1 1 1 0];
    % elseif numberOfOutlets==7
    % outletindex=[1 1 1 1 0 1 1 1];
    % elseif numberOfOutlets==8
    % outletindex=[1 1 1 1 1 1 1 1];
    % end

    % 如果使用遗传算法进行迭代，注释掉下面两行
    %load("GREI-outlet=1-42-4.232208.mat")
    %inp=population(6,:);

    BIM=dlmread('BIM.txt');  % 管道矩阵
    loops=dlmread('Loops.txt'); % 集水区矩阵
    condults=dlmread('[CONDUITS].txt'); % 管道上下游节点矩阵
    nodescordinates=dlmread('nodescordinates.txt');
    numberOfCatchmentArea = size(loops, 1);  % 18个集水区
    numberOfPipes = size(BIM, 1);            % 65个管道
    numberOfNodes= size(nodescordinates, 1); % 48个上下游节点
    cut=inp(2:numberOfCatchmentArea * 2 + 1);
    inppipe=inp(numberOfCatchmentArea * 2 + 2: numberOfCatchmentArea * 2 + 1 + size(outletindex, 2) + numberOfPipes); 

    for i=1:numberOfCatchmentArea      % 将管道能控制的集水区在矩阵中置为1
        for j=2:size(find(loops(i,:)),2)    
            BIM(loops(i,j),i+1)=1;          
        end
    end    

    BIM(:,numberOfCatchmentArea+2: numberOfCatchmentArea+3)=condults(:,2:3);   %将管道的上下节点进行替换
        
%% Allocate possible outlet  分配可能的出口

    % possibleRootNodes=[41 42 43 44 45 46 47 48];      
    % possibleRootpipes=[58 59 60 61 62 63 64 65];
    possibleRootNodes = firstRootNodes + (0:arrayLength-1); % 创建一个从n到n+length-1的数组 
    possibleRootpipes = firstRootpipes + (0:arrayLength-1); % 创建一个从n到n+length-1的数组 
    j = 1;
    for i=1:size(outletindex,2)
        if outletindex(1,i)==1
         Outlets(1,j) = possibleRootNodes(1,i);                       % 将出口记录下来
         j=j+1;
        else
        inppipe(1, size(outletindex,2) + possibleRootpipes(1,i)) = 0; % 将inppipe中不是出口对应的管道参数置为0
        end
    end

    Outlet = Outlets(1,1);  % 只取第一个出口
    dlmwrite('Outlets.txt',Outlets);

    for i=1:size(Outlets,2)
        num(i)=find(possibleRootNodes==Outlets(i)); % 出口的index
    end

    dlmwrite('Outpipes.txt',possibleRootpipes(num));% 将出口对应的管道记录下来
%% Cutting Algorithm  切割算法
        
    cdir=round(cut(1:numberOfCatchmentArea)); % Beta
    ncpipe=cut(numberOfCatchmentArea+1:2*numberOfCatchmentArea); % Alpha
    el=zeros(numberOfCatchmentArea,2);
    subCatchout=zeros(numberOfCatchmentArea,2);
    BIM2=BIM;

    for i=1:numberOfCatchmentArea % 每次循环代表一个集水区, i是集水区index，找到每个集水区对应切割管道
        numberOfPipesPerArea=sum(BIM(:,i+1)); % 每个集水区拥有管道数量
        cutIndex=round(ncpipe(i)*(numberOfPipesPerArea-1)+1); % This indicates that Copipe(th) in the column of Loop i is cut
        listOfPipesPerArea=find (BIM(:,1+i)>0); % 集水区i包含的管道组成的list
        allGone=isempty(listOfPipesPerArea);
        if allGone==1
            cutpipe(i)=0;
            break
        else
            cutpipe(i)=listOfPipesPerArea(cutIndex);
        end
        if cdir(i)==1
            upNode=BIM(cutpipe(i),numberOfCatchmentArea+3);
            downNode=BIM(cutpipe(i),numberOfCatchmentArea+2);
            BIM(cutpipe(i),numberOfCatchmentArea+3)=BIM(cutpipe(i),numberOfCatchmentArea+2);
        else
            upNode=BIM(cutpipe(i),numberOfCatchmentArea+2);
            downNode=BIM(cutpipe(i),numberOfCatchmentArea+3);
        end
        BIM(cutpipe(i),numberOfCatchmentArea+2) = numberOfNodes + i;
        % subcatchments and nodes for swmm
        subCatchout(i,1) = i + numberOfNodes + numberOfCatchmentArea;
        subCatchout(i,2) = numberOfNodes + i;
        % subCatchout(i,2)=Nup;Outlet of cathment
        nodescordinates(i+numberOfNodes,1)=i+numberOfNodes;
        nodescordinates(i+numberOfNodes,2)=(nodescordinates(upNode,2)+nodescordinates(downNode,2))*.5;
        nodescordinates(i+numberOfNodes,3)=(nodescordinates(upNode,3)+nodescordinates(downNode,3))*.5;

        el(i,1)=numberOfNodes+i;
        el(i,2)=BIM(cutpipe(i),numberOfCatchmentArea+4);
        BIM(cutpipe(i),2:numberOfCatchmentArea+1)=0; % When a pipe is selected from a loop to be cut, it is not included for other loops remained to be cut
        nodeOfEnd=BIM(cutpipe(i),numberOfCatchmentArea+3);
        BIM2(cutpipe(i),:)=0;
        q1=find (BIM2(:,numberOfCatchmentArea+3)==nodeOfEnd); % 下游节点是nodeOfEnd的管道
        q2=find (BIM2(:,numberOfCatchmentArea+2)==nodeOfEnd); % 上游节点是nodeOfEnd的管道
        while (length(q1)+length(q2))==1 % This guarantees that no node is ignored. Star problem
            if isempty(q1)
                NCutpipe(i)=q2;
                nodeOfEnd=BIM(NCutpipe(i),numberOfCatchmentArea+3);
            else
                NCutpipe(i)=q1;
                nodeOfEnd=BIM(NCutpipe(i),numberOfCatchmentArea+2);
            end
            BIM(NCutpipe(i),2:numberOfCatchmentArea+1)=0; % This guarantees that no node is ignored
            BIM2(NCutpipe(i),:)=0;
            q1=find (BIM2(:,numberOfCatchmentArea+3)==nodeOfEnd); % the pipes whose downstream node is Nend
            q2=find (BIM2(:,numberOfCatchmentArea+2)==nodeOfEnd);
        end
        %Check upstream of new node
    %     Nue(i)=Nup;
        BIM2(cutpipe(i),:)=0;
        q1=find (BIM2(:,numberOfCatchmentArea+3)==upNode); % the pipes whose downstream node is Nend
        q2=find (BIM2(:,numberOfCatchmentArea+2)==upNode); % the pipes whose upstream node is Nend
        while (length(q1)+length(q2))==1 % This guarantees that no node is ignored. Star problem
            if isempty(q1)
                NCutpipe(i)=q2;
                upNode=BIM(NCutpipe(i),numberOfCatchmentArea+3);
            else
                NCutpipe(i)=q1;
                upNode=BIM(NCutpipe(i),numberOfCatchmentArea+2);
            end
            BIM(NCutpipe(i),2:numberOfCatchmentArea+1)=0; % This guarantees that no node is ignored
            BIM2(NCutpipe(i),:)=0;
            q1=find (BIM2(:,numberOfCatchmentArea+3)==upNode); % the pipes whose downstream node is Nend
            q2=find (BIM2(:,numberOfCatchmentArea+2)==upNode);
        end
    end


        
%% Directing Graph
    Nnode= numberOfNodes + numberOfCatchmentArea; % Number of nodes in the spanning tree that is NP+1
    s = BIM(:,numberOfCatchmentArea+2);
    t = BIM(:,numberOfCatchmentArea+3);
    GG = graph(s,t);
    TR = shortestpathtree(GG,'all',Outlet);
    if size(TR.Edges,1)==size(BIM,1)
        %避免table2array(TR.Edges)排序
        BIM3(:,1:2)= table2array(TR.Edges);
        for i=1:numberOfPipes
            for j=1:numberOfPipes  
                if j < Outlet && BIM(i,numberOfCatchmentArea + 2) == BIM3(j,1) && BIM(i,numberOfCatchmentArea+3)~=BIM3(j,2) 
                    if BIM(i,numberOfCatchmentArea+3)< Outlet %需要替换所选出口的最小编号，只有一个出口时填这个出口，多个出口时填最小的编号
                        BIM(i,numberOfCatchmentArea+2:numberOfCatchmentArea+3)=BIM3(BIM(i,numberOfCatchmentArea+3),1:2);
                    else %一个出口时为'&&',多个出口时为'&"
                        BIM(i,numberOfCatchmentArea+2:numberOfCatchmentArea+3)=BIM3(BIM(i,numberOfCatchmentArea+3)-1,1:2);
                    end
                elseif j>=Outlet && BIM(i,numberOfCatchmentArea+2)==BIM3(j-1,1) && BIM(i,numberOfCatchmentArea+3)~=BIM3(j-1,2)
                    if BIM(i,numberOfCatchmentArea+3)< Outlet
                        BIM(i,numberOfCatchmentArea+2:numberOfCatchmentArea+3)=BIM3(BIM(i,numberOfCatchmentArea+3),1:2);
                    else
                        BIM(i,numberOfCatchmentArea+2:numberOfCatchmentArea+3)=BIM3(BIM(i,numberOfCatchmentArea+3)-1,1:2);
                    end
                end
            end
        end

        % BIM(:,numberOfCatchmentArea+2:numberOfCatchmentArea+3)= table2array(TR.Edges);
        s = BIM(:,numberOfCatchmentArea+2);
        t = BIM(:,numberOfCatchmentArea+3);
        sort=zeros(Nnode,1);
        INC=zeros(Nnode,Nnode);
        sort(1,1)=Outlet;
        Up=s(t==Outlet,1);
        k=size(Up,1);
        counter=1;
        for i=1:Nnode-1
            if k~=0
                for j=1:k
                    INC(i,j+counter)=1;
                    INC(counter+j,i)=1;
                    sort(counter+j,1)=Up(j,1);
                end
            end
            counter=counter+k;
            Down=sort(i+1,1);
            Up=s(t==Down,1);
            k=size(Up,1);
        end
        INC=[sort INC];
        pipes=[BIM(:,1) s t BIM(:,numberOfCatchmentArea+6) BIM(:,numberOfCatchmentArea+6)];
        p=BIM(:,1);
        EdgeTable = table([s t],p, ...
        'VariableNames',{'EndNodes' 'Pipes'});
        G = digraph(EdgeTable);
        plot(G,'EdgeLabel',G.Edges.Pipes)
        dlmwrite('subCatchout.txt',subCatchout);
        dlmwrite('[COORDINATES].txt',nodescordinates);
        dlmwrite('AhvazINC.txt',INC);
        dlmwrite('pipes.txt',pipes);
    else
        allGone=2;
    end
    if allGone==1
        obj=10e+12;
    elseif allGone==2
        obj=10e+12;
    else
        obj = Decentralize_hy(inppipe, dinp, flag);
    end
end

