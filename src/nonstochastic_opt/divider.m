function out = divider(NodeinnewPart,INC,Newrootp,Positions)%Update INC
Erorr=0;
SizenewInc=size(NodeinnewPart,2);
NewINC=zeros(SizenewInc,SizenewInc);
OldINC=INC;
NewPositions=zeros(SizenewInc,1);
NewPositions(1,1)=Newrootp;
newp=find(OldINC(Newrootp,:));
a=size(newp,2);
k=1;
sizetabu=1;
while k~=SizenewInc
for i=1:a
    NewPositions(sizetabu+i,1)=newp(1,i);
end
q=1;
for j=1:SizenewInc
    for i=1:a
    if NewPositions(j,1)==newp(1,i)
        temp(q,1)=j;
        q=q+1;
    end
    end
end
    
sizetabu=sizetabu+a;
for i=1:a
    NewINC(k,temp(i,1))=1;
    NewINC(temp(i,1),k)=1;   
end
k=k+1;
Upnodep=NewPositions(k,1);
if Upnodep==0
    Erorr=1;
    break
else
newp=find(OldINC(Upnodep,:));

t=0;
for i=1:size(newp,2)
    for j=1:SizenewInc
        if newp(1,i)==NewPositions(j,1)
            t=t+1;
            Tabu(1,t)=i;
        end
    end
end
 for i=1:t
     newp(:,Tabu(t))=[];
      Tabu=Tabu-1;
 end
a=size(newp,2);    
end
end
if Erorr==1
    out=zeros(size(NewINC,2),size(NewINC,2)+1);
    out(1, 1)=3000000;
else
out=zeros(size(NewINC,2),size(NewINC,2)+1);
out(:,1)=Positions(NewPositions);      
out(:,2:size(NewINC,2)+1)=NewINC; 
end
