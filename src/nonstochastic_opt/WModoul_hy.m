function SWMMInputW= WModoul_hy(pardesign) %
%%
    % pardesign(1,530)=.5;%Error in tabu!
    % pardesign(1,121)=1;%Keep at least one lid
    % pardesign(1,341:349)=0;% Outlet pipes
    %nproject=30;%project life time
    %inflatation=0.05; %%%%%%%%%%%%%%%%%%%%%%%%%！！！！！！！！！！！！！！！！！！！！通货膨胀率
    %discountrate=0.15;%%%%%%%%%%%%%%%%%%%%%%%%%！！！！！！！！！！！！！！！！！！！！折旧率
    %OMGray=.1;%Maintain and operation %%%%%%%%%！！！！！！！！！！！！！！！！！！！！灰色操作维护费率
%% MAX MIN Fix Layout
    Commercial=dlmread('Commercial.txt');%%%%%%%%%%%%%%%%%%%%%%%%%！！！！！！！！！！！！！！！！！！！！
    
    pipes=dlmread('pipes.txt');
    numberOfPipes = size(pipes, 1);
    for i = 1:numberOfPipes  % 将每个管道
        maxDiameters(i,1)=i;
        maxDiameters(i,2)=2;
    end
    MaxDiametersOptimal=dlmread('MaxDiametersOptimal.txt');
    for i=1:size(MaxDiametersOptimal,1)
        tf=find(MaxDiametersOptimal(i,1)==maxDiameters(:,1));
        if tf>0
            maxDiameters(MaxDiametersOptimal(i,1),2)=MaxDiametersOptimal(i,2);
        end
    end
    MaxDiametersOptimal=maxDiameters;
    MinDiametersOptimal(1:numberOfPipes,2)=Commercial(1,1);

% ZonesLID=dlmread('Zones.txt'); 

%% [SUBCATCHMENTS]
    %   pardesign=ones(1,530);
    %Hmax= 5.1255  Relia= 99.9 Cost= 9.1858e+04
    input_args=load('INCFinal.mat'); %这个是保存的数据矩阵文件
    input_args=input_args.INCFinal;
    Roughness=0.01;
    Subs=dlmread('Subcatchout.txt');
    Initial=dlmread('[SUBCATCHMENTS].txt');
    NodeCoo=dlmread('[COORDINATES].txt');
    NN=size(NodeCoo,1);%Num. Nodes
    myformat = '%1s%1s%16s%17s%17s%10s%10s%9s%7s%s%[^\n\r]';
    InputfileOrg = regexp( fileread('AhvazNull.inp'), '\n', 'split');
    Inputfile=InputfileOrg;
    ToSearch=regexp( fileread('[SUBCATCHMENTS] - T.txt'), '\n', 'split');
    tf=strcmp(Inputfile,ToSearch);%对应位置元素进行比较，相同为1，否则为0
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
    fid = fopen('LIMatF2.inp', 'w');  %这个.inp文件是什么？
    fprintf(fid, '%s\n', Inputfile{:});
    fclose(fid);

%%   Pipe Sizing
    input=input_args;
    NumberofParts=size(input,2);
    INCPipes=zeros(size(pipes,1),size(pipes,1)+1);
    INCPipes(:,1)=pipes(:,1);
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
             k=1;
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
    for i=1:size(MaxDiametersOptimal,1)
    Dmax(i)=min(max(Commercial(:,1)),MaxDiametersOptimal(i,2));
    end
    for i=1:size(MaxDiametersOptimal,1)
    Dmin(i)=max(DminMarket,MinDiametersOptimal(i,2));
     end
    %inp=[0,0.390888442539489,0.103955079430394,0.622993080958774,0.503627884556948,0.471453403912240,0.835240001801894,0.313076960670168,0.187073336870826,0.594860191768767,0.115407834042954,0.802932757750101,0.764389475330716,0.935371751364088,0.490442253953495,0.0728713824102303,0.278744101935837,0.463667455687445,0.917795310269657,0.542536742655418,0.947986504890022,0.0438208374642465,0.925252158345484,0.886922613781590,0.579758759243151,0.246020271698743,0.732340374428769,0.257537637193149,0.742458391525520,0.203220802799784,0.366209238175715,0.491045476354124,0.995099444989281,0.787960157972270,0.130621442902985,0.430269875720873,0.789292340120335,0.0777601254359405,0.604478639139869,0.0256341168795431,0.849508766971832,0.744230526696204,0.474516062199755,0.858942300965342,0.151601779702189,0.982931178008908,0.188430681258928,0.780149177876810,0.101581706491908,0.178130373383259,0.620995173146943,0.174343184582008,0.467717393260654,0.118533470574614,0.605755208135679,0.955307919579994,0.962208627483814,0.404185033785811,0.465708950072574,0.870052003538002,0.132959608596945,0.928500411825342,0.623496315003533,0.294649331868839,0.540941559262017,0.220423545834707,0.833343879300162,0.999061641196606,0.436863536967204,0.369777516490067,0.748829285187346,0.708355587139842,0.381725654934090,0.795600857758253,0.673277690463935,0.436707212177091,0.710826027048161,0.670078478648967,0.379588875575907,0.406318475384054,0.471985886045457,0.759682025935570,0.00254830765629865,0.575341993574652,0.531899632500369,0.436771100230556,0.730206872802076,0.678311998486003,0.588313972467694,0.664567275175626,0.282182216266248,0.622408708194465,0.182774950591139,0.883301088885522,0.459712217686848,0.305974209816001,0.220390038945197,0.374384312510660,0.908531575168861,0.431298189675353,0.921878398022710,0.440667672744379,0.0764387868852643,0.179548846589030,0.155670673949150,0.944676700306408,0.772197921006145,0.856667839222629,0.843387827597576,0.489079641028223,0.989318735235628,0.510858763164238,0.529199117195722,0.144012629530434,0.987977371398504,0.695936662040377,0.475794406022022,0.814993616334701,0.960406628099114,0.128042546031227,0.922462145263208,0.700344225302851,0.578219752021250,0.745379943895285,0.0990184360676850,0.143601870666306,0.694368514613798,0.859914863656914,0.649533579666797,0.410158736523979,0.854381983463181,0.403248587982267,0.561230212134243,0.719781472480003,0.393478283971185,0.888116239240273,0.929118911451127,0.0188893272318917,0.433498544188483,0.817204177721758,0.833637920373324,0.233976528770582,0.531151735316362,0.481070858060488,0.470814818278173,0.739073724610374,0.907216325296704,0.478978461982829,0.150861794347969,0.862729064418875,0.846807468700126,0.842367522113452,0.623646533727110]
    %inppipe=inp(60:153);
    %pardesign=inppipe(7:size(inppipe,2));
    input=pardesign;
    %%单独运行的话需要复制变量
    INCPipesTemp=INCPipes;
    INCPipesTemp(:,1)=[];
    for i=1:NumberofParts
        TempInfopipes=Infopipes{i};
        Npipes=size(TempInfopipes,1);
        k=1;
        for j=Npipes:-1:1
            PipeName=TempInfopipes(j,1);
            Uppipes=find(INCPipesTemp(PipeName,:));
            s=size(Uppipes,2);
            if s==0
                D=input(1,PipeName)*(Dmax(PipeName)-Dmin(PipeName))+Dmin(PipeName);%D=input(1,k)*(DmaxMarket-DminMarket)+DminMarket;

                for f=1:size(Commercial,1)-1
                    if D>=Commercial(f,1) &&  D<(Commercial(f+1,1)+Commercial(f,1))*0.5
                        D=Commercial(f,1);
                        Smin=Commercial(f,2);%minimum slope from the list
                        break
                    elseif  D>=(Commercial(f+1,1)+Commercial(f,1))*0.5 &&  D<=Commercial(f+1,1)
                            D=Commercial(f+1,1);
                            Smin=Commercial(f,2);%minimum slope from the list
                            break
                    else
                    end
                end
                TempInfopipes(j,5)=D;
                TempInfopipes(j,13)=Smin;
            else 
                %Dmin=Dmin(PipeName);
                for z=1:s
                    up=find(TempInfopipes(:,1)==Uppipes(1,z));
                    Dup=TempInfopipes(up,5);
                    DminTel=max(Dmin(PipeName),Dup);
                    Dmin(PipeName)=DminTel;
                end
                D=input(1,PipeName)*(Dmax(PipeName)-DminTel)+DminTel;%D=input(1,k)*(DmaxMarket-Dmin)+Dmin;

                for f=1:size(Commercial,1)-1
                    if D>=Commercial(f,1) &&  D<(Commercial(f+1,1)+Commercial(f,1))*0.5
                       D=Commercial(f,1);
                       Smin=Commercial(f,2);%minimum slope from the list
                       break
                    elseif D>=(Commercial(f+1,1)+Commercial(f,1))*0.5 &&  D<=Commercial(f+1,1)
                           D=Commercial(f+1,1);
                           Smin=Commercial(f,2);%minimum slope from the list
                           break
                    else
                    end
                end
                TempInfopipes(j,5)=D;
                TempInfopipes(j,13)=Smin;
            end     
            k=k+1;%variable countr
        end
        Infopipes{i}=TempInfopipes;
    end

    % - Slopes Node elevations and pipe ofset
    Cover=0.8;
    NodeEle=dlmread('NodeEle.txt');
    outlets=zeros(NumberofParts,2);
    for i=1:NumberofParts
        TempInfopipes=Infopipes{i};
        TempInconodes=input_args{i};
        Npipes=size(TempInfopipes,1);
        Nnodes=size(TempInconodes,1);
        k=1;
        for j=Nnodes:-1:2
            NodeName=TempInconodes(j,1);
            Conectednodes=find(TempInconodes(j,:));
            Conectednodes(1)=[];
            Conectednodes=Conectednodes-1;
            sizeconected=size(Conectednodes,2);
            Downnode=Conectednodes(1,1);
            DownnodeName=TempInconodes(Downnode);
            upnodes=Conectednodes(1,2:sizeconected);
            upnodesName=TempInconodes(upnodes,1)';
            sizeupnodes=size(upnodes,2);
            if sizeupnodes~=0
               uppipes=zeros(1,sizeupnodes);
               for l=1:sizeupnodes
                   uppipes(1,l)=Nodes2pipe(upnodesName(1,l),NodeName,pipes);
               end
            else
            end
            pipeinhand=Nodes2pipe(NodeName,DownnodeName,pipes);
            pipeinhandLoc=find(TempInfopipes(:,1)==pipeinhand);
            if sizeupnodes==0
                InEU=NodeEle(NodeName,2)-Cover-TempInfopipes(pipeinhandLoc,5);
                InED=min(InEU-TempInfopipes(pipeinhandLoc,13)*TempInfopipes(pipeinhandLoc,4),NodeEle(DownnodeName,2)-1);
                TempInfopipes(pipeinhandLoc,6)=InEU;
                TempInfopipes(pipeinhandLoc,7)=InED;
            else
                UppipesDInv=zeros(1,sizeupnodes);
                UppipesDCrown=zeros(1,sizeupnodes);
                for k=1:sizeupnodes
                    UppipeLo(1,k)=find(TempInfopipes(:,1)==uppipes(1,k));
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
    % uselesspipes=dlmread('uselesspipes.txt');
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

    %断1根120个散点；断2-118根各200个散点；断119根120个散点；断120根1个散点
    %除掉废管

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
        Inputfile{1,LineToEdit+i} = sprintf(myformat,Temp1,'CIRCULAR     ',Temp2);%这个代码和;LIDModoul里的不一样
    end

    fid = fopen('LIMatF2.inp', 'w'); %清空上面的数据
    fprintf(fid, '%s\n', Inputfile{:});%重新写入
    fclose(fid);
    %% Cost calculation
    NP=size(CONDUITS,1);
     for i=1:NP
        CONDUITS(i,9)=18-CONDUITS(i,6);%H manholes
        CONDUITS(i,10)=18-.5*(CONDUITS(i,6)+CONDUITS(i,7));%dhm miangin omghe kargozaari har looleh
     end
    dH=CONDUITS(:,10);
    dHm=CONDUITS(:,9);
    L=CONDUITS(:,4);
    maxalow=12;
    for i=1:NP
        D(i)=CONDUITS(i,5);
        if D(i)==.2
                CostP(i)=(3.4881*dH(i) - 2.2851)*L(i);
                CostM(i)=(29.952*dHm(i) + 54.912);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
        elseif D(i)==.25
                CostP(i)=(3.5112*dH(i) - 2.1031)*L(i);
                CostM(i)=(33.696*dHm(i) + 58.656);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
        elseif D(i)==.3
                CostP(i)=(3.3674*dH(i) - 1.4712)*L(i);
                CostM(i)=(38.359*dHm(i) + 58.11);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
        elseif D(i)==.35
                CostP(i)= (3.3773*dH(i) - 1.3078)*L(i);
                CostM(i)=(42.027*dHm(i) + 62.212);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
        elseif D(i)==.4
                CostP(i)=(3.5669*dH(i) - 0.907)*L(i);
                CostM(i)=(45.663*dHm(i) + 66.456);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
        elseif D(i)==.53
                CostP(i)=(3.6239*dH(i) + 0.0551)*L(i);
                CostM(i)=(52.845*dHm(i) + 75.374);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
        elseif D(i)==.6
                CostP(i)=(3.8548*dH(i) + 1.5862)*L(i);
                CostM(i)=(56.39*dHm(i) + 80.047);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
        elseif D(i)==.8
                CostP(i)= (4.3059*dH(i) + 3.8577)*L(i);
                CostM(i)=(59.904*dHm(i) + 84.864);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
         elseif D(i)==1
                CostP(i)= (4.65*dH(i) + 8.12)*L(i);
                CostM(i)=(73.65*dHm(i) + 105.56);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
         elseif D(i)==1.20
                CostP(i)= (5.11*dH(i) + 10.84)*L(i);
                CostM(i)=(79.4*dHm(i) + 113.6);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
         elseif D(i)==1.50
                CostP(i)= (5.73*dH(i) + 15.55)*L(i);
                CostM(i)=(91.2*dHm(i) + 129.86);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
         elseif D(i)==2.00
                CostP(i)= (6.78*dH(i) + 23.38)*L(i);
                CostM(i)=(110.904*dHm(i) + 159.32);
                if dH(i)>maxalow
                     CostP(i)=CostP(i)*dH(i);
                end
         elseif D(i)==2.40
                CostP(i)= (7.7*dH(i) + 28.38)*L(i);
                CostM(i)=(130.904*dHm(i) + 188.32);
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

    CashflowGray=zeros(1,nproject);
    CostGray=(sum(CostP(1:size(CONDUITS,1)))+sum(CostM(1:size(CONDUITS,1)))); % cost function
    for i=1:nproject
        CashflowGray(1,i)=OMGray*CostGray*(1/(1+discountrate)^i);
    end
    TotalOMGray=sum(CashflowGray);

    Scost=(CostGray+TotalOMGray)*1000000/42105*2; %9880774.95+CostGray+TotalOMGray

    %   SWMMInputW=1000*H;
%% RUN SWMM and read data
    ! swmm5.exe  LIMatF2.inp  LIMatF2.rpt 
    filename = 'LIMatF2.rpt';
    delimiter = {'..','\t',',',' ',';'};
    startRow = 58;%
    endRow = 58;%
    formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false); %%读入大文件
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID); 
    LIMatF = [dataArray{1:end-1}];
    a=LIMatF(1,1);
    if  a=='Flooding'
        Flood=str2num(LIMatF(1,4));
    else 
        Flood=999;    
    end
        Vflooding=Flood;

    if Vflooding<0.001
        SWMMInputW=Scost;
    else 
        SWMMInputW=Scost*(1+10000*Vflooding);%SWMMInputW=Scost*(1+100*Vflooding)
    end    




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



