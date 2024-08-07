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
