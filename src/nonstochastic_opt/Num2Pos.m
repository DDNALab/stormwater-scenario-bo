function out = Num2Pos(inp,NumofNodes,Positions)%Convert position of nodes to real number of them
sizeinp=size(inp,2);
out=zeros(1,sizeinp);
for i=1:sizeinp
for j=1:NumofNodes
    if inp(1,i)==Positions(j,1);
        out(1,i)=j;
    end
end
end