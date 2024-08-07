function out = pathfinder(NumofNodes,INC,newrootp,oldrootp)%Convert position of nodes to real number of them
% PathP=zeros(1,NumofNodes);
Pathp(1,1)=newrootp;
UpNodep=newrootp;
k=1;
% if pathp(1,1)==0
%     out=0;
% else
while UpNodep~=oldrootp
    DownNodeP=UpNodep;
    for i=1:DownNodeP
    if INC(DownNodeP,i)==1
        k=k+1;
        UpNodep=i;
        Pathp(1,k)=UpNodep;
    break
    end
    end

end

out=Pathp;