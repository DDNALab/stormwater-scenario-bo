function SWMMLID= LIDModoul(pardesign,dinp)%灰色时为SWMMInputW= WModoul(pardesign)
%This modoule optimize LID and Gray network
% StructuralResilience=StructuralResilience*-.01;
%%
% pardesign(1,530)=.5;%Error in tabu!
% dinp=randi([0,1],1,88*2+29*2+29*2);%Keep at least one lid %前119个是管径 dinp(1,119*2+1)
% pardesign(1,341:349)=0;% Outlet pipes
%% MAX MIN Fix Layout
Commercial=dlmread('Commercial.txt');%%%%%%%%%%%%%%%%%%%%%%%%%！！！！！！！！！！！！！！！！！！！！
pipes=dlmread('pipes.txt');
loops=dlmread('Loops.txt'); % 集水区矩阵
numberOfCatchmentArea = size(loops, 1);  % 18个集水区
numberOfPipes = size(pipes, 1);
MaxDiameters=zeros(1,numberOfPipes);  %替换管道数量
for i=1:numberOfPipes %替换管道数量
    MaxDiameters(i,1)=i;
    MaxDiameters(i,2)=2;
end
MaxDiametersOptimal=dlmread('MaxDiametersOptimal.txt');
for i=1:size(MaxDiametersOptimal,1)
    tf=find(MaxDiametersOptimal(i,1)==MaxDiameters(:,1));
    if tf>0
        MaxDiameters(MaxDiametersOptimal(i,1),2)=MaxDiametersOptimal(i,2);
    else
    end
end
MaxDiametersOptimal=MaxDiameters;
% Min Diameter n smaller than current (Max)
% NDimMin=2;
% for i=1:530
%     row=find(Commercial(:,1)==MaxDiametersOptimal(i,2));
%     MinDiametersOptimal(i,1)=i;
%     MinDiametersOptimal(i,2)=Commercial(max(1,row-NDimMin),1);
% end

 MinDiametersOptimal(1:numberOfPipes,2)=Commercial(1,1); %替换管道数量

% ZonesLID=dlmread('Zones.txt'); 
%% Initialization 初始化
%[SUBCATCHMENTS]
input_args=load('INCFinal.mat'); 
input_args=input_args.INCFinal;
Roughness=0.01;
Subs=dlmread('Subcatchout.txt');
Initial=dlmread('[SUBCATCHMENTS].txt');
NodeCoo=dlmread('[COORDINATES].txt');
NN=size(NodeCoo,1);%Num. Nodes
myformat = '%1s%1s%16s%17s%17s%10s%10s%9s%7s%s%[^\n\r]';
InputfileOrg = regexp( fileread('AhvazNull.inp'), '\n', 'split');    %在AhvazNull.inp文件中可以更改LID类型和参数，包括surface，soil，storage,drain等
 
Inputfile=InputfileOrg;
  
ToSearch=regexp( fileread('[SUBCATCHMENTS] - T.txt'), '\n', 'split');
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;
for i=1:size(Subs,1)
    Initial(i,1)=Subs(i,1);
    Initial(i,3)=Subs(i,2);
    Temp=num2str(Initial(i,:));
    Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp);
end

% [SUBAREAS]         OUTLET 
SUBAREAS=dlmread('[SUBAREAS].txt');
formatSpec ='%1s%1s%19s%10s%12s%11s%9s%15s%s%[^\n\r]';
ToSearch=regexp( fileread('[SUBAREAS]-T.txt'), '\n', 'split');
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;
for i=1:size(Subs,1)
    SUBAREAS(i,1)=Subs(i,1);
    Temp=num2str(SUBAREAS(i,:));
    Temp2='    OUTLET';
    Inputfile{1,LineToEdit+i} = sprintf(formatSpec,Temp,Temp2);
end
% [INFILTRATION]
INFILTRATION=dlmread('[INFILTRATION].txt');
myformat = '%1s%1s%18s%11s%12s%s%[^\n\r]';
ToSearch=regexp( fileread('[INFILTRATION]-T.txt'), '\n', 'split');
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;
for i=1:size(Subs,1)
    INFILTRATION(i,1)=Subs(i,1);
    Temp=num2str(INFILTRATION(i,:));
    Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp);
end
% 
% [Polygons]
Polygons=dlmread('[Polygons].txt');
myformat = '%1s%1s%23s%19s%s%[^\n\r]';
ToSearch=regexp( fileread('[Polygons]-T.txt'), '\n', 'split');
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;
for i=1:size(Subs,1)
    Polygons(i,1)=Subs(i,1);
    Temp=num2str(Polygons(i,:));
    Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp);
end
%fid = fopen('LIMatF2.inp', 'w'); 用于打开文件
%fprintf(fid, '%s\n', Inputfile{:});
%fclose(fid);

%%   Pipe Sizing
pipes=dlmread('pipes.txt');
input=input_args;
NumberofParts=size(input,2);
INCPipes=zeros(size(pipes,1),size(pipes,1)+1);
INCPipes(:,1)=pipes(:,1);
Infopipes=cell(1,NumberofParts);
for i=1:NumberofParts% This generate the information about pipes
    nodes=input{i};
    nodesnumber=nodes(:,1);
    nodesinc=nodes(:,2:size(nodes,2));
    Numberofpipes=size(nodesinc,1)-1;
    INCpipestemp=zeros(Numberofpipes,12);%INCpipestemp=zeros(Numberofpipes,Numberofpipes+1);
     for j=Numberofpipes:-1:1
         upnode=nodesnumber(j+1,1);
         conections=find(nodesinc(j+1,:));
         nconection=size(conections,2)-1;
         downnode=nodesnumber(conections(1,1),1);
         pipenumber=Nodes2pipe(upnode,downnode,pipes);
         INCpipestemp(j,1)=pipenumber;
         INCpipestemp(j,2)=upnode;
         INCpipestemp(j,3)=downnode;
         INCpipestemp(j,4)=pipes(pipenumber,5);
%          k=1;
         if nconection~=0
             for k=1:nconection
                 conectednode=nodesnumber(conections(1,k+1),1);
                 upstreampipe=Nodes2pipe(upnode,conectednode,pipes);
                 INCPipes(pipenumber,upstreampipe+1)=1;%incidence matrix of pipes (Upstram pipes of each pipe= 1 in row
             end
         else
         end
     end
     
     Infopipes{i}=INCpipestemp;%%%%%info pipes
     %PipeNumber   UpNode   DoNode   L   D   S   LoUpIn   LoDoIn   GUpIn GDoIn   Ofset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% Pipe sizing Dmin Dmax alfa
% 1- Diameters
DminMarket=min(Commercial(:,1));
Dmax=zeros(1,size(MaxDiametersOptimal,1));
Dmin=zeros(1,size(MaxDiametersOptimal,1));
for i=1:size(MaxDiametersOptimal,1)
    Dmax(i)=min(max(Commercial(:,1)),MaxDiametersOptimal(i,2));
end
for i=1:size(MaxDiametersOptimal,1)
    Dmin(i)=max(DminMarket,MinDiametersOptimal(i,2));
end
%inp=[0.1637293381,0.7525,0.89609,0.36038,0.84615,0.32562,0.62271,0.089745,0.38655,0.35678,0.96883,0.091189,0.97793,0.99278,0.016465,0.15321,0.90748,0.8508,0.40737,0.519,0.058923,0.49504,0.67715,0.1315,0.096616,0.89587,0.56158,0.68053,0.43642,0.87085,0.13302,0.18268,0.29546,0.54452,0.22319,0.3868,0.62047,0.95917,0.26251,0.68872,0.66354,0.4935,0.7102,0.12943,0.8714,0.16352,0.29291,0.165,0.28456,0.31222,0.99462,0.072861,0.52153,0.46454,0.091129,0.28566,0.79071,0.377,0.76117,0.71443,0.053512,0.005515,0.74096,0.89709,0.62696,0.54526,0.86321,0.14576,0.94489,0.24388,0.87156,0.7483,0.15978,0.80692,0.50213,0.43047,0.016617,0.39632,0.2684,0.62209,0.59429,0.19214,0.64744,0.056617,0.26363,0.60102,0.11847,0.79303,0.80001,0.74747,0.44478,0.45013,0.89861,0.53563,0.17418,0.44209,0.12216,0.17175,0.55307,0.70881,0.93957,0.62592,0.79436,0.95126,0.82438,0.36286,0.40544,0.74481,0.38771,0.97994,0.63536,0.70657,0.67013,0.74998,0.74894,0.89979,0.52187,0.92712,0.90628,0.53619,0.64608,0.032703,0.88309,0.064107,0.87516,0.99722,0.62553,0.44169,0.081576,0.38833,0.70424,0.34375,0.8246,0.7351,0.34733,0.52496,0.76033,0.28201,0.93234,0.30436,0.94772,0.36033,0.31214,0.91473,0.95294,0.089491,0.71216,0.32706,0.93478,0.86127,0.70404,0.14112,0.70549]
%inppipe=inp(60:153);
%pardesign=inppipe(7:size(inppipe,2));%单独运行需补上这三行代码
input=pardesign;
INCPipesTemp=INCPipes;
INCPipesTemp(:,1)=[];
%load('dinp.mat')
%dinp(:,177)=1
%dinp=dinp(81,:)
%for i=1:204
%    dinp=dinp(i,:)
   
for i=1:NumberofParts
    TempInfopipes=Infopipes{i};
    Npipes=size(TempInfopipes,1);
    k=1;
    for j=Npipes:-1:1 %遍历每一根管
        PipeName=TempInfopipes(j,1);
        Uppipes=find(INCPipesTemp(PipeName,:));%找出所有上管管名
        s=size(Uppipes,2);
        if s==0 %无上管的管
            D=input(1,PipeName)*(Dmax(PipeName)-Dmin(PipeName))+Dmin(PipeName);%分配每根管的D
            for f=1:size(Commercial,1)-1
                if D>=Commercial(f,1) &&  D<(Commercial(f+1,1)+Commercial(f,1))*0.5
                    D=Commercial(f,1);
%                     Smin=Commercial(f,2);%minimum slope from the list
                    break
                elseif  D>=(Commercial(f+1,1)+Commercial(f,1))*0.5 &&  D<=Commercial(f+1,1)
                    D=Commercial(f+1,1);
%                     Smin=Commercial(f+1,2);%minimum slope from the list
                    break
                else
                end
            end
%dinp=randi([0,1],204,292)
   %%%这里选最优解dinp
           %%%%得到最优解集并找到最优解后，在这里更改dinp得到韧性和成本
%加入dinp判断是否减管径      pipename  1-130  
           q=find(Commercial(:,1)==D);
           if q==2 && dinp(1,PipeName*2)==0 && dinp(1,PipeName*2-1)==1 %10减一次
               D=Commercial(q-1,1);
           elseif q==3 && dinp(1,PipeName*2)==0 && dinp(1,PipeName*2-1)==1 %10减一次
               D=Commercial(q-1,1);
           elseif q==3 && dinp(1,PipeName*2)==1 && dinp(1,PipeName*2-1)==0 %01减二次
               D=Commercial(q-2,1);
           elseif q>=4 && dinp(1,PipeName*2)==0 && dinp(1,PipeName*2-1)==1 %10减一次
               D=Commercial(q-1,1);
           elseif q>=4 && dinp(1,PipeName*2)==1 && dinp(1,PipeName*2-1)==0 %01减二次
               D=Commercial(q-2,1);
           elseif q>=4 && dinp(1,PipeName*2)==0 && dinp(1,PipeName*2-1)==0 %00减三次
               D=Commercial(q-3,1);
           end
           Smin=Commercial(Commercial(:,1)==D,2);%找出该管的最小坡度

                TempInfopipes(j,5)=D;
                TempInfopipes(j,13)=Smin;
                
        else %有上管的管
            
%             %限制最小管径一定≥上管管径
            for z=1:s
                up=TempInfopipes(:,1)==Uppipes(1,z);%寻找上管在TempInfopipes中的行数
                Dup=TempInfopipes(up,5);%获取上管管径
                DminTel=max(Dmin(PipeName),Dup);%最小管径一定≥上管管径
                Dmin(PipeName)=DminTel;%获得最小管径
            end
            D=input(1,PipeName)*(Dmax(PipeName)-DminTel)+DminTel;%得到灰色最优基础上的管径后，①开启减管径，②减完检查是否满足≥上管管径
            %D=inp中该管数值*（Simpleobj限制的最大管径-刚刚获得的最小管径）+刚刚获得的最小管径；例如管109D=0.7446*（0.6-0.25）+0.25=0.5106
%            
            %D四舍五入并获取相应的最小坡度Smin
            for f=1:size(Commercial,1)-1 
                if D>=Commercial(f,1) &&  D<(Commercial(f+1,1)+Commercial(f,1))*0.5
                    D=Commercial(f,1);
%                     Smin=Commercial(f,2);%minimum slope from the list
                    break
                elseif  D>=(Commercial(f+1,1)+Commercial(f,1))*0.5 &&  D<=Commercial(f+1,1)
                    D=Commercial(f+1,1);
%                     Smin=Commercial(f+1,2);%minimum slope from the list
                    break
                else
                end
           end %得D=0.53
           
           %加入dinp判断是否减管径
           D2=D;
           q=find(Commercial(:,1)==D);
           if q==2 && dinp(1,PipeName*2)==0 && dinp(1,PipeName*2-1)==1 %10减一次
               D2=Commercial(q-1,1);
           elseif q==3 && dinp(1,PipeName*2)==0 && dinp(1,PipeName*2-1)==1 %10减一次
               D2=Commercial(q-1,1);
           elseif q==3 && dinp(1,PipeName*2)==1 && dinp(1,PipeName*2-1)==0 %01减二次
               D2=Commercial(q-2,1);
           elseif q>=4 && dinp(1,PipeName*2)==0 && dinp(1,PipeName*2-1)==1 %10减一次
               D2=Commercial(q-1,1);
           elseif q>=4 && dinp(1,PipeName*2)==1 && dinp(1,PipeName*2-1)==0 %01减二次
               D2=Commercial(q-2,1);
           elseif q>=4 && dinp(1,PipeName*2)==0 && dinp(1,PipeName*2-1)==0 %00减三次
               D2=Commercial(q-3,1);
           end
           
           %再次限制最小管径一定≥上管管径
           Dshang=zeros(1,s);
            for z=1:s
                up=TempInfopipes(:,1)==Uppipes(1,z);%寻找上管在TempInfopipes中的行数
                Dup=TempInfopipes(up,5);%获取上管管径
                DminTel=max(Dmin(PipeName),Dup);%选择最大的上管管径（Dmin(PipeName)-其他上管管径；Dup-该上管管径）
                Dmin(PipeName)=DminTel;%更新上管管径
                Dshang(z)=DminTel;
            end
            if D2>=max(Dshang) %如果减后满足≥上管管径
                D=D2; %管径变为减后管径D2
            end
           Smin=Commercial(Commercial(:,1)==D,2);%找出该管的最小坡度
%             D=input(1,PipeName)*(Dmax(PipeName)-DminTel)+DminTel;%分配每根管的D
            %D=inp中该管数值*（Simpleobj限制的最大管径-刚刚获得的最小管径）+刚刚获得的最小管径；例如管109D=0.7446*（0.6-0.25）+0.25=0.5106
%            Smin=Commercial(find(Commercial(:,1)==TempInfopipes(j,5)),2);%找出该管的最小坡度
           TempInfopipes(j,5)=D;
           TempInfopipes(j,13)=Smin;
        end     
        k=k+1;%variable countr
    end
    %main的for循环针对每一根管减一次管径，判断是否flooding，若flooding则保留，不flooding则减
%     if TempInfopipes(decreasepipe,5)>0.25
%         w=find(Commercial(:,1)==TempInfopipes(decreasepipe,5));
%         TempInfopipes(decreasepipe,5)=Commercial(w-1,1);
%     end    
    Infopipes{i}=TempInfopipes;
end
 


% - Slopes Node elevations and pipe ofset
Cover=0.8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!
NodeEle=dlmread('NodeEle.txt');
outlets=zeros(NumberofParts,2);
for i=1:NumberofParts
    TempInfopipes=Infopipes{i};
    TempInconodes=input_args{i};
    %Npipes=size(TempInfopipes,1);
    Nnodes=size(TempInconodes,1);
    %k=1;
    for j=Nnodes:-1:2
        NodeName=TempInconodes(j,1);
        Conectednodes=find(TempInconodes(j,:));
        Conectednodes(1)=[];
        Conectednodes=Conectednodes-1;
        sizeconected=size(Conectednodes,2);
        Downnode=Conectednodes(1,1);
        DownnodeName=TempInconodes(Downnode);%找出此节点的下节点
        upnodes=Conectednodes(1,2:sizeconected);%找出此节点的上节点
        upnodesName=TempInconodes(upnodes,1)';
        sizeupnodes=size(upnodes,2);
        if sizeupnodes~=0
            uppipes=zeros(1,sizeupnodes);
            for l=1:sizeupnodes
                uppipes(1,l)=Nodes2pipe(upnodesName(1,l),NodeName,pipes);
            end
        else
        end
        pipeinhand=Nodes2pipe(NodeName,DownnodeName,pipes);%找出节点下的管
        pipeinhandLoc=find(TempInfopipes(:,1)==pipeinhand);%找出管在TempInfopipes中的位置
        if sizeupnodes==0
            InEU=NodeEle(NodeName,2)-Cover-TempInfopipes(pipeinhandLoc,5);%18-0.8-那根管径
            InED=min(InEU-TempInfopipes(pipeinhandLoc,13)*TempInfopipes(pipeinhandLoc,4),NodeEle(DownnodeName,2)-1);%min（InEU-管最小坡度*长度，初始高度-1）
            TempInfopipes(pipeinhandLoc,6)=InEU;
            TempInfopipes(pipeinhandLoc,7)=InED;
        else
            UppipesDInv=zeros(1,sizeupnodes);
            UppipesDCrown=zeros(1,sizeupnodes);
            for k=1:sizeupnodes
                UppipeLo(1,k)=find(TempInfopipes(:,1)==uppipes(1,k));%57
                UppipesDInv(1,k)=TempInfopipes(UppipeLo(1,k),7);
                UppipesDCrown(1,k)=UppipesDInv(1,k)+TempInfopipes(UppipeLo(1,k),5);
            end
            MinCrown=min(UppipesDCrown);
            InEU=MinCrown-TempInfopipes(pipeinhandLoc,5);
            InED=min(InEU-TempInfopipes(pipeinhandLoc,13)*TempInfopipes(pipeinhandLoc,4),NodeEle(DownnodeName,2)-1);
            TempInfopipes(pipeinhandLoc,6)=InEU;
            TempInfopipes(pipeinhandLoc,7)=InED;            
            for k=1:sizeupnodes
                 ofset=TempInfopipes(UppipeLo(1,k),7)-InEU;
                 if ofset<0.001
                     ofset=0;
                 end
                 TempInfopipes(UppipeLo(1,k),8)=ofset;
                 TempInfopipes(UppipeLo(1,k),7)=InEU;
             end
        end
    end
    outlets(i,1)=TempInfopipes(1,3);
    outlets(i,2)=TempInfopipes(1,7);
    TempInfopipes(:,13)=[];
    Infopipes{i}=TempInfopipes;
end

%除掉废管
uselesspipes=dlmread('uselesspipes.txt');
% for i=1:size(Infopipes,2)
%     for k=1:size(uselesspipes,1)
%         for u=1:size(Infopipes{i},1)
%             if uselesspipes(k,1)==Infopipes{i}(u,1)
%                 Infopipes{i}(u,:)=0;
%             end
%         end
%     end
%     Infopipes{i}(find(Infopipes{i}(:,1)==0),:)=[];
% end

% Writing nodes information
Junctions=zeros(NN-NumberofParts,6);
k=1;
for i=1:NumberofParts
    Temp=Infopipes{i};
    nnodes=size(Temp,1);
    Junctions(k:k+nnodes-1,1)=Temp(:,2);
    Junctions(k:k+nnodes-1,2)=Temp(:,6);
    k=k+nnodes;
end
for j=1:k-1
Junctions(j,3)=NodeEle(Junctions(j,1),2)-Junctions(j,2);
end
Junctions=sortrows(Junctions,1);
% Junctions(find(Junctions(:,1)==0),:)=[];

%% Nodes
%[COORDINATES]
COORDINATES=dlmread('[COORDINATES].txt');
myformat = '%s%s%s%[^\n\r]';
ToSearch=regexp( fileread('[COORDINATES]-T.txt'), '\n', 'split');
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;
for i=1:NN
    Temp=num2str(COORDINATES(i,:));
    Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp);
end

%[JUNCTIONS]
ToSearch=regexp( fileread('[JUNCTIONS]-T.txt'), '\n', 'split');
myformat = '%2s%16s%11s%11s%11s%11s%s%[^\n\r]';
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;
for i=1:size(Junctions,1)
    Temp=num2str(Junctions(i,:));
    Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp);
end
% [OUTFALLS]
% ;;Name           Elevation  Type       Stage Data       Gated    Route To        
% ;;-------------- ---------- ---------- ---------------- -------- ----------------
% 56               0          FREE                        NO  
ToSearch=regexp( fileread('[OUTFALLS]-T.txt'), '\n', 'split');
myformat = '%s%s%s%s%[^\n\r]';
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;

for i=1:NumberofParts
    Temp=num2str([outlets(i,1) outlets(i,2)]);
    Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp,'          FREE                        NO');
end

%% Pipes
% [CONDUITS]
% ;;Name           From Node        To Node          Length     Roughness  InOffset   OutOffset  InitFlow   MaxFlow   
% ;;-------------- ---------------- ---------------- ---------- ---------- ---------- ---------- ---------- ----------
NumofPipes=0;
for i=1:NumberofParts
NumofPipes=NumofPipes+size(Infopipes{i},1);
end
CONDUITS=zeros(NumofPipes,12);
k=1;
for i=1:NumberofParts
    Temp=Infopipes{i};
    npipes=size(Temp,1);
    CONDUITS(k:k+npipes-1,:)=Temp;
    k=k+npipes;
end
ToSearch=regexp( fileread('[CONDUITS]-T.txt'), '\n', 'split');
myformat = '%2s%17s%17s%18s%12s%8s%11s%11s%11s%s%[^\n\r]';
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;
LineToEdit2=LineToEdit;
Line2=Line;
Inputfile2=Inputfile;
CONDUITS2=CONDUITS;
myformat2=myformat;
Initial2=Initial;
for i=1:size(CONDUITS,1)
    Temp=[CONDUITS(i,1)	CONDUITS(i,2)	CONDUITS(i,3)	CONDUITS(i,4)	Roughness	CONDUITS(i,9)	CONDUITS(i,8)	CONDUITS(i,10)	CONDUITS(i,11)];
    Temp=num2str(Temp);
    Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp);
end

% [XSECTIONS]
% ;;Link           Shape        Geom1            Geom2      Geom3      Geom4      Barrels    Culvert   
% ;;-------------- ------------ ---------------- ---------- ---------- ---------- ---------- ----------
myformat = '%2s%23s%6s%17s%11s%11s%11s%s%[^\n\r]';
ToSearch=regexp( fileread('[XSECTIONS]-T.txt'), '\n', 'split');
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+2;
for i=1:size(CONDUITS,1)
    Temp1=CONDUITS(i,1);
    Temp1=num2str(Temp1);
    Temp2=[CONDUITS(i,5)	0	0	0	1];
    Temp2=num2str(Temp2);
    Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp1,'CIRCULAR     ',Temp2);
end
% fid = fopen('LIMatF2.inp', 'w');
% fprintf(fid, '%s\n', Inputfile{:});
% fclose(fid);

% [LID_USAGE]
% ;;Subcatchment   LID Process      Number  Area       Width      InitSat    FromImp    ToPerv     RptFile                  DrainTo         
% ;;-------------- ---------------- ------- ---------- ---------- ---------- ---------- ---------- ------------------------ ----------------
%LID1 IT AREA=20 With=2 FromImp=45  LID2 RB AREA=1.5  FromImp=35
%717    1	15	20	2	0	45	0 
InitSat=0;
AreaBC=200; %每个BC单元200㎡   %这个可能要换，一个单元200平对永庆坊来说可能太大了 大概50㎡
WidthBC=10;
AreaPP=150; %每个PP单元150㎡   %这个大概37.5㎡
WidthPP=10;
ToPervBC=1; %所有溢流返回至透水区域
ToPervPP=1; %所有溢流返回至透水区域
LIDProcessBC=1;
LIDProcessPP=2;

%dinp(1,257:301)控制FromImpRB和FromImpBC
% fromImpBC=100*dinp(1,249:291); %处理不渗透面积比
% fromImpPP=100-fromImpBC; %处理不渗透面积比

%dinp(1,167:256)控制NumofBCperHec和NumofPPperHec，约束：BC+PP面积≤汇水区面积10%
% NumofBCperHec=5*dinp(1,163:205);%珠江新城最大83800*10%ha的汇水区面积
% NumofPPperHec=20/3*dinp(1,206:248);%珠江新城最大83800*10%ha的汇水区面积
% for i=1:43  %替换集水区数量
%     if dinp(1,size(CONDUITS,1)*2+43+i*2-1)==1 && dinp(1,size(CONDUITS,1)*2+43+i*2)==1 %11乘1/4
%         NumofBCperHec(i)=1/4*5;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*2-1)==0 && dinp(1,size(CONDUITS,1)*2+43+i*2)==1 %01乘2/4
%         NumofBCperHec(i)=2/4*5;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*2-1)==1 && dinp(1,size(CONDUITS,1)*2+43+i*2)==0 %10乘3/4
%         NumofBCperHec(i)=3/4*5;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*2-1)==0 && dinp(1,size(CONDUITS,1)*2+43+i*2)==0 %00乘4/4
%         NumofBCperHec(i)=5;
%     end
% end
% for i=1:43
%     if dinp(1,size(CONDUITS,1)*2+43+i*3-2)==0 && dinp(1,size(CONDUITS,1)*2+43+i*3-1)==0 && dinp(1,size(CONDUITS,1)*2+43+i*3)==0 %11乘1/4
%         NumofPPperHec(i)=0/7*20/3;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*3-2)==0 && dinp(1,size(CONDUITS,1)*2+43+i*3-1)==0 && dinp(1,size(CONDUITS,1)*2+43+i*3)==1%11乘1/4
%         NumofPPperHec(i)=1/7*20/3;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*3-2)==0 && dinp(1,size(CONDUITS,1)*2+43+i*3-1)==1 && dinp(1,size(CONDUITS,1)*2+43+i*3)==0%11乘1/4
%         NumofPPperHec(i)=2/7*20/3;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*3-2)==1 && dinp(1,size(CONDUITS,1)*2+43+i*3-1)==0 && dinp(1,size(CONDUITS,1)*2+43+i*3)==0%11乘1/4
%         NumofPPperHec(i)=3/7*20/3;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*3-2)==0 && dinp(1,size(CONDUITS,1)*2+43+i*3-1)==1 && dinp(1,size(CONDUITS,1)*2+43+i*3)==1%11乘1/4
%         NumofPPperHec(i)=4/7*20/3;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*3-2)==1 && dinp(1,size(CONDUITS,1)*2+43+i*3-1)==1 && dinp(1,size(CONDUITS,1)*2+43+i*3)==0%11乘1/4
%         NumofPPperHec(i)=5/7*20/3;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*3-2)==1 && dinp(1,size(CONDUITS,1)*2+43+i*3-1)==0 && dinp(1,size(CONDUITS,1)*2+43+i*3)==1%11乘1/4
%         NumofPPperHec(i)=6/7*20/3;
%     elseif dinp(1,size(CONDUITS,1)*2+43+i*3-2)==1 && dinp(1,size(CONDUITS,1)*2+43+i*3-1)==1 && dinp(1,size(CONDUITS,1)*2+43+i*3)==1%11乘1/4
%         NumofPPperHec(i)=7/7*20/3;
%     end
% end

% CostRB=19.423;%M Rials 每个RB多少钱 
% CostIT=113.4;%M Rials 每单元IT多少钱  
% OMRB=0;%Operation Maintain
% OMIT=.04;%Operation Maintain

% NumofBC=0;
% NumofPP=0;
myformat = '%2s%23s%6s%17s%11s%11s%11s%s%[^\n\r]';  %%指定输出格式，%s是输出字符串
ToSearch=regexp( fileread('[LID_USAGE]-T.txt'), '\n', 'split');  % regexp(     'slipe')用于在指定位置将字符串分开
tf=strcmp(Inputfile,ToSearch);
Line=find(tf);
LineToEdit=Line+3;
CapBC=zeros(1,size(Initial,1));
CapPP=zeros(1,size(Initial,1));
AreaofBC=zeros(1,size(Initial,1));
AreaofPP=zeros(1,size(Initial,1));

%Without Zones  这一部分就开始判断集水区中LID的面积
for i=1:size(Initial,1)
%      返回矩阵initial的行数也就是集水区的数量
%将initial替换成除建筑外的面积
% if dinp(i+119*2)==1 %1有LID
%      if pardesign(i+size(pipes,1))>.5
     Subcatchment=Initial(i,1);%%%集水区编号  90——   从junction（54）+outfall（6）+集水区（29）开始编号
     Area=Initial(i,4);%汇水区面积
     %Usable area = Area-Floorarea%(这个是需要导入数据的)
%          imperv=Initial(i,5);%不渗透比
     %Areaimperv=imperv*0.01*Area*10000;%不渗透部分面积
%          Areaperv=(1-imperv*0.01)*Area*10000;%渗透部分面积
    %%我们的目的：将除建筑面积外的剩余面积（Area-Floor area)全部改为PP或者BC，且PP>BC。
    %面积占(Area-Floor area)的百分之多少    在这里加入Area-Floor area数据
    %编码就是00=0, 01=
%%%%替换管道数量    i=18  &&&从131列开始
    if dinp(1,numberOfPipes*2+i*2)==0 && dinp(1,numberOfPipes*2+i*2-1)==0 %00 0%  %%%dinp(1,119*2+i*2)----
       LIDSizeBC=0;%%%判断第178，131列-166,233列
    elseif dinp(1,numberOfPipes*2+i*2)==1 && dinp(1,numberOfPipes*2+i*2-1)==0 %01 3.33%
       LIDSizeBC=Area*5/3;%1.6667
    elseif dinp(1,numberOfPipes*2+i*2)==0 && dinp(1,numberOfPipes*2+i*2-1)==1 %10 6.67%
       LIDSizeBC=Area*5*2/3;%3.3333
    elseif dinp(1,numberOfPipes*2+i*2)==1 && dinp(1,numberOfPipes*2+i*2-1)==1 %11 10%
       LIDSizeBC=Area*5;
    end%%%%判断第236列-292列
    if dinp(1,numberOfPipes*2+numberOfCatchmentArea+numberOfCatchmentArea+i*2)==0 && dinp(1,numberOfPipes*2+numberOfCatchmentArea+numberOfCatchmentArea+i*2-1)==0 %00 0%
       LIDSizePP=0;
    elseif dinp(1,numberOfPipes*2+numberOfCatchmentArea+numberOfCatchmentArea+i*2)==1 && dinp(1,numberOfPipes*2+numberOfCatchmentArea+numberOfCatchmentArea+i*2-1)==0 %01 3.33%
       LIDSizePP=Area*20/3/3;%2.2222
    elseif dinp(1,numberOfPipes*2+numberOfCatchmentArea+numberOfCatchmentArea+i*2)==0 && dinp(1,numberOfPipes*2+numberOfCatchmentArea+numberOfCatchmentArea+i*2-1)==1 %10 6.67%
       LIDSizePP=Area*20/3*2/3;%4.4444
    elseif dinp(1,numberOfPipes*2+numberOfCatchmentArea+numberOfCatchmentArea+i*2)==1 && dinp(1,numberOfPipes*2+numberOfCatchmentArea+numberOfCatchmentArea+i*2-1)==1 %11 10%
       LIDSizePP=Area*20/3;%6.6667
    end %%PP的计算在BC 的计算上×了4/3  原代码是小于集水区面积的10%，那永庆坊就是小于（Area-floor area)的百分百
    if LIDSizeBC~=0 || LIDSizePP~=0  %~=不等于  这里额外加入判断条件？PP>BC？
         nBC=round(LIDSizeBC); %nIT=round(Area*NumofITperHec) %子集水区面积*每公顷子集水区的BC比 取整
         nPP=round(LIDSizePP); %nRB=round(Area*NumofRBperHec) %子集水区面积*每公顷子集水区的PP比
         %BC+PP面积不能大于该汇水区的10%
         BCbili=AreaBC*nBC/(AreaBC*nBC+AreaPP*nPP); %BC在总LID中所占的比例
         PPbili=AreaPP*nPP/(AreaBC*nBC+AreaPP*nPP);
         if nBC*AreaBC+nPP*AreaPP>Area*1000 %如果加起来大于10%，按10%乘以BCPP比例  nBC*AreaBC+nPP*AreaPP>(Area- Floor AREA)*10000
             nBC=round((Area*1000*BCbili)/AreaBC);  %nBC=round((Usable area*10000*BCbili)/AreaBC)
             nPP=round((Area*1000*PPbili)/AreaPP);   %%nPP=round((Usable area*10000*PPbili)/AreaPP)
             BCbili=AreaBC*nBC/(AreaBC*nBC+AreaPP*nPP);   %这里的AreaBC*nBC+AreaPP*nPP会变成1000以下， BC\PPbili不变，但NBC和NPP变
             PPbili=AreaPP*nPP/(AreaBC*nBC+AreaPP*nPP);
         end 
    %          NumofBC=nBC+NumofBC; %NumofRB=nRB+NumofRB
    %          NumofPP=nPP+NumofPP; %NumofIT=NumofIT+nIT
         FromImpBC=BCbili*100;%FromImp按两者面积比
         FromImpPP=PPbili*100;
    %          FromImpBC=fromImpBC(i);
    %          FromImpPP=fromImpPP(i);
         if nBC==0
             FromImpBC=0;
             FromImpPP=100;
         end
         if nPP==0
             FromImpPP=0;
             FromImpBC=100;
         end
         Temp1=[Subcatchment	LIDProcessBC	nBC	AreaBC	WidthBC InitSat FromImpBC ToPervBC]; %Temp1=[Subcatchment	LIDProcessIT	nIT	AreaIT	WidthIT InitSat FromImpIT ToPervIT]; 
         Temp2=[Subcatchment	LIDProcessPP	nPP	AreaPP	WidthPP InitSat FromImpPP ToPervPP]; %Temp2=[Subcatchment	LIDProcessRB	nRB	AreaRB	WidthRB InitSat FromImpRB ToPervRB];
         Temp1=num2str(Temp1);
         Temp2=num2str(Temp2);
         Inputfile{1,LineToEdit} = sprintf(myformat,Temp1);
         Inputfile{1,LineToEdit+1} = sprintf(myformat,Temp2);
         LineToEdit=LineToEdit+2;
         %下面就是算每一个集水区LID的成本了
         AreaofBC(1,i)=nBC*AreaBC; %列出每个子汇水区的BC面积（公顷）
         CapBC(1,i)=75.6*AreaofBC(1,i)+15*(AreaofBC(1,i)^0.5); %126938.47
         AreaofPP(1,i)=nPP*AreaPP;
         CapPP(1,i)=34.6*AreaofPP(1,i)+15*(AreaofPP(1,i)^0.5); %126938.47   %这里计算成本  
         %^矩阵幂
%      else
     end
end%%540行的for in 循环

fid = fopen('LIMatF2.inp', 'w');%打开文件清空内容
fprintf(fid, '%s\n', Inputfile{:});  %将数据写入文本文件
fclose(fid);
% 50 years Storm %1                INTENSITY 0:01     1.0      TIMESERIES 20
LineToEdit=51;
% sprintf函数将数据格式化为字符串或字符向量
Inputfile{1,LineToEdit} = sprintf(myformat,'1                VOLUME 0:01     1.0      TIMESERIES 5yChicago-6h');  %%设计降雨修改
fid = fopen('LIMatF50.inp', 'w');
fprintf(fid, '%s\n', Inputfile{:});
fclose(fid);

%% Cost calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%广州灰色价格
NP=size(CONDUITS,1);
 for i=1:NP
    CONDUITS(i,9)=18-CONDUITS(i,6);%H manholes
    CONDUITS(i,10)=18-.5*(CONDUITS(i,6)+CONDUITS(i,7));%dhm miangin omghe kargozaari har looleh
 end
dH=CONDUITS(:,10);
dHm=CONDUITS(:,9);
L=CONDUITS(:,4);
maxalow=12;
CostP=zeros(1,NP);
CostM=zeros(1,NP);

for i=1:NP
    D(i)=CONDUITS(i,5);
    if D(i)==.2
        CostP(i)=(3.4881*dH(i) - 2.2851)*L(i);
        CostM(i)=29.952*dHm(i) + 54.912;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==.25
        CostP(i)=(3.5112*dH(i) - 2.1031)*L(i);
        CostM(i)=33.696*dHm(i) + 58.656;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==.3
        CostP(i)=(3.3674*dH(i) - 1.4712)*L(i);
        CostM(i)=38.359*dHm(i) + 58.11;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==.35
        CostP(i)= (3.3773*dH(i) - 1.3078)*L(i);
        CostM(i)=42.027*dHm(i) + 62.212;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==.4
        CostP(i)=(3.5669*dH(i) - 0.907)*L(i);
        CostM(i)=45.663*dHm(i) + 66.456;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==.53
        CostP(i)=(3.6239*dH(i) + 0.0551)*L(i);
        CostM(i)=52.845*dHm(i) + 75.374;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==.6
        CostP(i)=(3.8548*dH(i) + 1.5862)*L(i);
        CostM(i)=56.39*dHm(i) + 80.047;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==.8
        CostP(i)= (4.3059*dH(i) + 3.8577)*L(i);
        CostM(i)=59.904*dHm(i) + 84.864;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==1
        CostP(i)= (4.65*dH(i) + 8.12)*L(i);
        CostM(i)=73.65*dHm(i) + 105.56;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==1.20
        CostP(i)= (5.11*dH(i) + 10.84)*L(i);
        CostM(i)=79.4*dHm(i) + 113.6;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==1.50
        CostP(i)= (5.73*dH(i) + 15.55)*L(i);
        CostM(i)=91.2*dHm(i) + 129.86;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==2.00
        CostP(i)= (6.78*dH(i) + 23.38)*L(i);
        CostM(i)=110.904*dHm(i) + 159.32;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==2.40
        CostP(i)= (7.7*dH(i) + 28.38)*L(i);
        CostM(i)=130.904*dHm(i) + 188.32;
            if dH(i)>maxalow
                 CostP(i)=CostP(i)*dH(i);
            end
    elseif D(i)==0.01
        CostP(i)= 0;
        CostM(i)=0;

    end
end

nproject=30;%project life time
% inflatation=0.05; %通货膨胀率
discountrate=0.02;%折旧率
OMGray=.1;%Maintain and operation %灰色操作维护费率
OMPP=0.04;
OMBC=0.08;

CashflowGray=zeros(1,nproject);
CashflowPP=zeros(1,nproject);
CashflowBC=zeros(1,nproject);
CostGray=(sum(CostP(1:NP))+sum(CostM(1:NP)))*1000000/42105*2; % cost function
for i=1:nproject
    CashflowGray(1,i)=OMGray*CostGray*(1/(1+discountrate)^i);
    CashflowPP(1,i)=OMPP*sum(CapPP)*(1/(1+discountrate)^i);
    CashflowBC(1,i)=OMBC*sum(CapBC)*(1/(1+discountrate)^i);
end
TotalOMGray=sum(CashflowGray);
TotalOMPP=sum(CashflowPP);
TotalOMBC=sum(CashflowBC);
% TotalLCCBC=sum(LCCBC);%7122517.54
CostPP=sum(CapPP);%2219041.96
CostBC=sum(CapBC);%2219041.96
Scost=CostPP+TotalOMPP+CostBC+TotalOMBC+CostGray+TotalOMGray; %9880774.95+CostGray+TotalOMGray
% Scost=TotalCostRB+TotalCostIT+CostGray+TotalOMIT+TotalOMGray;
%533.33㎡，75不透水率，不透水面积为400，5%不透水面积为20平，IT一单元20平，IT一单元+30年+RB一个=133.4+231.7654+19.423=8659.7刀
%BC=245.85n+51*根号n，PP=76.12n+33*根号n，n=533.33*5%=26.66666平，BC+PP=9019.7刀
dlmwrite('COSTPP.txt',CostPP)  %% 这里计算可以写入
dlmwrite('COSTBC.txt',CostBC)
dlmwrite('COSTTotalOMPP.txt',TotalOMPP)
dlmwrite('COSTTotalOMBC.txt',TotalOMBC)
dlmwrite('COSTCostGray.txt',CostGray)
dlmwrite('COSTTotalOMGray.txt',TotalOMGray)
dlmwrite('Scost.txt',Scost)
% Scost=Scost/10;%Rial to Toman
%% Reliability T=2 
% Reliability 2y
! swmm5.exe  LIMatF2.inp  LIMatF2.rpt 
filename = 'LIMatF2.rpt';
delimiter = {'..','\t',',',' ',';'};
startRow = 40;   %%这个值是什么意思？
endRow = 59;                                                                                   
formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
LIMatF = [dataArray{1:end-1}];  %%这个矩阵的问题
Floodingline=find(LIMatF(:,1)=='Flooding');      
Totalrunoffline=find(LIMatF(:,1)=='Surface');
Vflooding=str2double(LIMatF(Floodingline,4));    % 第18行，第4列，vflooding=22.1030    str2double将字符串转换为双精度值   Floodingline=18  与集水区坡度，不透水面比例有
Totalrunoff5y=str2double(LIMatF(Totalrunoffline,4));
 
%下面这个判断无法进行或者说判断为SWMMLID=[10000000000000 1
%10000000000000]，这样就无法进行OPR和TBR的计算，进而导致遗传算法中值不变
if Vflooding<=0.001    %%%0.001      
    % Reliability 2y
    ! swmm5.exe  LIMatF50.inp  LIMatF50.rpt 
    filename = 'LIMatF50.rpt';
    delimiter = {'..','\t',',',' ',';'};
    startRow = 40;
    endRow = 59;
    formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);  %%%从文本文件或字符串读取格式化数据。
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    LIMatF = [dataArray{1:end-1}];
    Floodingline=find(LIMatF(:,1)=='Flooding');
    Totalrunoffline=find(LIMatF(:,1)=='Surface');
    Vflooding50=str2double(LIMatF(Floodingline,4));
    Totalrunoff50y=str2double(LIMatF(Totalrunoffline,4));
    Res_TBR=1-(Vflooding50/161.247/1.0182);
    Res_OPR=OPRA(LineToEdit2,Line2,Inputfile2,CONDUITS2,myformat2,Initial2,dinp,uselesspipes);
   %%%如果加入BC随时间而衰退的内容，那直接相应的缩减BC面积比例可行吗？？
    Res=-100000000*((Res_TBR*Res_OPR)^0.5);   
    dlmwrite('ResTBR.txt',Res_TBR)
    dlmwrite('ResOPR.txt',Res_OPR)
    Totalrunoff=10000000*(Totalrunoff5y*Totalrunoff50y)^0.5;
    SWMMLID=[Scost Res Totalrunoff];   
    %dlmwrite ('multiobj.txt',CostPP,CostBC,TotalOMPP,TotalOMBC,CostGray,TotalOMGray,Scost,Res_OPR,Res_TBR);
else 
    SWMMLID=[10000000000000 1 10000000000000];
end    
end


% ! swmm5.exe  LIMatF2.inp  LIMatF2.rpt 
% filename = 'LIMatF2.rpt';
% delimiter = {'..','\t',',',' ',';'};
% startRow = 58;
% endRow = 58;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF = [dataArray{1:end-1}];
% a=LIMatF(1,1);
% if a=='Flooding'
% Flood=str2num(LIMatF(1,4));
% else 
% Flood=999;    
% end
% Vflooding=Flood;
% 
% startRow = 44;
% endRow = 44;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF2 = [dataArray{1:end-1}];
% b=LIMatF2(1,1);
% if b=='Surface' 
% Runoff=str2num(LIMatF2(1,4));
% else 
% Runoff=999;    
% end
% Vrunoff=Runoff;
% Reliability2=1-(Vflooding/Vrunoff);
% 
% if Vflooding<=0.001
% Reliability2=1-(Vflooding/Vrunoff);
% else 
% Reliability2=0;
% end
% 
% if Reliability2==0
%     SWMMLID=Scost*(1+100*Vflooding)
% else
%% Reliability T=10
% Reliability 10y
% ! swmm5.exe  LIMatF10.inp  LIMatF10.rpt 
% filename = 'LIMatF10.rpt';
% delimiter = {'..','\t',',',' ',';'};
% startRow = 58;
% endRow = 58;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF = [dataArray{1:end-1}];
% a=LIMatF(1,1);
% if a=='Flooding'
% Flood=str2num(LIMatF(1,3));
% else 
% Flood=999;    
% end
% Vflooding=Flood;
% 
% startRow = 44;
% endRow = 44;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF2 = [dataArray{1:end-1}];
% b=LIMatF2(1,1);
% if b=='Surface' 
% Runoff=str2num(LIMatF2(1,4));
% else 
% Runoff=999;    
% end
% Vrunoff=Runoff;
% 
% % if Vflooding<=Vrunoff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%是否还应判断10年一遇最低地表径流量？
% Reliability10=1-(Vflooding/Vrunoff);
% % else 
% % Reliability10=0;
% % end
% 
% % if Reliability10==0
% %     SWMMLID=[Scost*(1+100*Vflooding) 0]
% % else
% %% Reliability 20y
% % 
% ! swmm5.exe  LIMatF20.inp  LIMatF20.rpt 
% filename = 'LIMatF20.rpt';
% delimiter = {'..','\t',',',' ',';'};
% startRow = 58;
% endRow = 58;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF = [dataArray{1:end-1}];
% a=LIMatF(1,1);
% if a=='Flooding'
% Flood=str2num(LIMatF(1,3));
% else 
% Flood=999;    
% end
% Vflooding=Flood;
% 
% startRow = 44;
% endRow = 44;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF2 = [dataArray{1:end-1}];
% b=LIMatF2(1,1);
% if b=='Surface' 
% Runoff=str2num(LIMatF2(1,4));
% else 
% Runoff=999;    
% end
% Vrunoff=Runoff;
% Reliability20=1-(Vflooding/Vrunoff);
% 
% % Vflooding=Flood/90;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%90是总径流量
% % Reliability20=1-Vflooding;
% 
% %% Resilience 25y
% ! swmm5.exe  LIMatF25.inp  LIMatF25.rpt 
% filename = 'LIMatF25.rpt';
% delimiter = {'..','\t',',',' ',';'};
% startRow = 58;
% endRow = 58;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF = [dataArray{1:end-1}];
% a=LIMatF(1,1);
% if a=='Flooding'
% Flood=str2num(LIMatF(1,3));
% else 
% Flood=999;    
% end
% Vflooding=Flood;
% 
% startRow = 44;
% endRow = 44;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF2 = [dataArray{1:end-1}];
% b=LIMatF2(1,1);
% if b=='Surface' 
% Runoff=str2num(LIMatF2(1,4));
% else 
% Runoff=999;    
% end
% Vrunoff=Runoff;
% Reliability25=1-(Vflooding/Vrunoff);
% % 
% % Vflooding=Flood/170;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%170是总径流量
% % Reliability25=1-Vflooding;
% %% Resilience T=50 and Structures
% ! swmm5.exe  LIMatF50.inp  LIMatF50.rpt 
% filename = 'LIMatF50.rpt';
% delimiter = {'..','\t',',',' ',';'};
% startRow = 58;
% endRow = 58;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF = [dataArray{1:end-1}];
% a=LIMatF(1,1);
% if a=='Flooding'
% Flood=str2num(LIMatF(1,3));
% else 
% Flood=999;    
% end
% Vflooding=Flood;
% 
% startRow = 44;
% endRow = 44;
% formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% LIMatF2 = [dataArray{1:end-1}];
% b=LIMatF2(1,1);
% if b=='Surface' 
% Runoff=str2num(LIMatF2(1,4));
% else 
% Runoff=999;    
% end
% Vrunoff=Runoff;
% Reliability50=1-(Vflooding/Vrunoff);


% Vflooding=Flood/170;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%170是总径流量
% Reliability50=1-Vflooding;
%% Environmental Sustainability 
% EnvironmentalSustainability=(NumofRB+NumofIT)/(5698+7589);
% EnvironmentalSustainability=(NumofBC+NumofPP)/(267+356);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5698+7589是最大LID数量

%% Objectives
% TotalSustainability=-10000000*(EnvironmentalSustainability*Reliability50*StructuralResilience*Reliability25*Reliability20*Reliability10*Reliability2)^(1/6); % StructuralResilience灰色弹性、每种暴雨情景权重
% SWMMLID=Scost;

% end
% end
%
  function Accordingpipe = Nodes2pipe(upnode,downnode,pipes)
%   pipes=dlmread('pipes.txt');
  X=pipes(:,2);
  Y=pipes(:,3);
  Accordingpipe=find(X==upnode & Y==downnode);
  TF = isempty(Accordingpipe);
  if TF==1
      Accordingpipe=find(X==downnode & Y==upnode);
  else
  end
  end
