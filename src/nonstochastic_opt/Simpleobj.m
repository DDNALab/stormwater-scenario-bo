function cost = Simpleobj(INCFinal)
%Simple objective for costs and resilience
%   Calculate area connected to each pipe

Subcatchout=dlmread('Subcatchout.txt');
SubcatchmentArea=dlmread('SubcatchmentArea.txt');
pipes=dlmread('pipes.txt');
pipes(:,5)=0;%pipes(:,5) Area directly connected to each upstream pipe
for i=1:size(Subcatchout,1)
    for j=1:size(pipes,1)
        if Subcatchout(i,2)==pipes(j,2)
            pipes(j,5)=SubcatchmentArea(i,2);
        end
    end
end


itt=size(INCFinal,2);
% pipes(:,6)=pipes(:,4);
pipes(:,6)=pipes(:,5);
for i=1:itt
    numnodes=size(INCFinal{i},1);
    INC=INCFinal{i};
    nodes=INC(:,1);
    INC(:,1)=[];
     for j=numnodes:-1:2
        Currentnode=nodes(j,1);
        conections=find(INC(j,:));
        downnode=nodes(conections(1,1));
        if size(conections,2)~=1
        conections(1)=[];
        for k=1:size(conections,2)
            upnodes(1,k)=nodes(conections(1,k),1);
        end
        else
            upnodes=0;
        end
        
        %Currentpipe
        a=[Currentnode downnode];
        for z=1:size(pipes,1)
            b(1,1)=pipes(z,2);
            b(1,2)=pipes(z,3);
            c(1,1)=pipes(z,3);
            c(1,2)=pipes(z,2);
            if a==b | a==c
                Currentpipe=z;
                break
            end
        end
        %UpPipes
        if upnodes~=0
            Uppipe=zeros(size(upnodes,2),1);
        for y=1:size(upnodes,2)
            a=[Currentnode upnodes(1,y)];
            for z=1:size(pipes,1)
            b(1,1)=pipes(z,2);
            b(1,2)=pipes(z,3);
            c(1,1)=pipes(z,3);
            c(1,2)=pipes(z,2);
            if a==b | a==c
                Uppipe(y,1)=z;
                break
            end
            end
        end
        else
            Uppipe=0;
        end

        %QT
        
        g=1;
     
        while Uppipe(g,1)~=0 
            pipes(Currentpipe,6)=pipes(Uppipe(g,1),6)+pipes(Currentpipe,6);
            g=g+1;
             if g>size(Uppipe,1)
                 break
             end
             
        end
            
        
     end
     
    
end
%Simple Cost
c=pipes(:,4).*(pipes(:,6).^0.5);
 cost(1,1)=sum(c);
 pipes(:,7)=2;
 %Area Diammeter 为了减少运算量先算出小管 区分不同的汇水区面积和不透水率
% uselesspipes=zeros(size(pipes,1),1);
 for i=1:size(pipes,1)
     if pipes(i,6)==0
         pipes(i,7)=0.25;
%          uselesspipes(i,1)=pipes(i,1)
%      elseif pipes(i,6)<=1
%          pipes(i,7)=0.4;
     elseif pipes(i,6)<=1.5
         pipes(i,7)=0.8;
     elseif pipes(i,6)<=5
         pipes(i,7)=1;
     elseif pipes(i,6)<=8
         pipes(i,7)=1.2;
     elseif pipes(i,6)<=13
         pipes(i,7)=1.5;
     else
         pipes(i,7)=2;
     end
 end
%增加无用管矩阵
% uselesspipes(find(uselesspipes==0))=[];
% dlmwrite('uselesspipes.txt',uselesspipes);
dlmwrite('MaxDiametersOptimal.txt',[pipes(:,1) pipes(:,7)]);
% DC
%  DC=size(INCFinal,2);
% if DC>0
%  cost(1,3)=DC;
% else
%     cost(1,3)=100;
% end

% % Structure Resilience
% weights=[1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1]';
% Npipes=size(pipes,1);
% TotalArea=sum(pipes(:,5));
% categories=zeros(1,10);
% pipes(:,7)=pipes(:,6)/TotalArea;
% for j=1:10
% for i=1:Npipes %#ok<ALIGN>
%     if pipes(i,7)>((j-1)/10)&&pipes(i,7)<=(j/10)
%         categories(1,j)=1+categories(1,j);
%     end
%     end
% end
% cost(1,2)=-100*(categories(1,1)/Npipes)*categories(1,2:10)*weights(2:10,1)/sum(categories(1,2:10));
% Structure Resilience 2
Npipes=size(pipes,1);
TotalArea=sum(pipes(:,5));
pipes(:,7)=1-(pipes(:,6)/TotalArea);
NRes90=Npipes;
summation=0;
for i=1:Npipes
    if pipes(i,7)<0.9 
        NRes90=NRes90-1;
        summation=summation+pipes(i,7);
    end
end
Resilience=-100*(NRes90/Npipes)*(summation/(Npipes-NRes90));
cost(1,2)=Resilience;

