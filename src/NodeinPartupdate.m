function out = NodeinPartupdate(Itteration,NodeinPart,contribiutingpart,NodeinnewPart)%update node in part matrix
%z=2时Itteration=z=2，NodeinPart=3x17 double，contribiutingpart=1,NodeinnewPart=[4,7,11,16]
for i=1:size(NodeinnewPart,2)             %1                  %2                  %3                   %4
    temp=NodeinnewPart(1,i);              %4                  %7                  %11                  %16
    NodeinPart(contribiutingpart,temp)=0; %NodeinPart(1,4)=0  %NodeinPart(1,7)=0  %NodeinPart(1,11)=0  %NodeinPart(1,16)=0
    NodeinPart(Itteration+1,temp)=1;      %NodeinPart(3,4)=1  %NodeinPart(3,7)=1  %NodeinPart(3,11)=1  %NodeinPart(3,16)=1
end
out =NodeinPart; 
%[1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,0,0; %oldoldpart
% 0,0,0,0,0,0,0,1,1,1,0,1,1,1,0,0,1; %oldpart
% 0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0] %newpart
%问题：z=2时NodeinPart还是3×17矩阵，但应该是3×10矩阵