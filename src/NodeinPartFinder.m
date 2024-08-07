function NodeinnewPart = NodeinPartFinder(INC,Cutpipeupnodep,NumofNodes)%Provide NodeinPart function
savenode=Cutpipeupnodep;    
NewInc=INC;
for i=1:NumofNodes%make under main diammeter zero
    for j=1:i
       NewInc(i,j)=0;
    end
end
 
NewNodesPosition=find(NewInc(Cutpipeupnodep,:));
AddedNodes=NewNodesPosition;

NumofNew=size(NewNodesPosition,2);
Sum=NumofNew;
b=1;
while b~=0
 num=0;
for i=1:NumofNew
Temp=find(NewInc(AddedNodes(1,i),:));
N=size(Temp,2);
NewNodesPosition(1,Sum+1:Sum+N)=Temp;
AddedNodestemp(1,num+1:num+N)=Temp;
num=num+N;
Sum=Sum+N;
end
AddedNodes=AddedNodestemp;
AddedNodestemp=0;
NumofNew=num;
b=NumofNew;
end
a=size(NewNodesPosition,2); 
b=NewNodesPosition;
NewNodesPosition(1,1)=savenode;
NewNodesPosition(2:a+1)=b;



NodeinnewPart=NewNodesPosition;
