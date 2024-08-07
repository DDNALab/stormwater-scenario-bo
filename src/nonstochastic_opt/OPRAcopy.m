function Res_OPR=OPRA(LineToEdit2,Line2,Inputfile,CONDUITS,myformat,Initial,dinp,uselesspipes)
%% 鏂?1鏍?120涓暎鐐癸紱鏂?2-118鏍瑰悇202涓暎鐐癸紱鏂?119鏍?120涓暎鐐癸紱鏂?120鏍?1涓暎鐐广?傛瘡绉嶇粍鍚堥兘瑕佽繍琛屼竴娆WMM锛屾瘡涓猧閮借SWMM202娆?
%鈶犳寜D鎺掑簭鈫扗paixupipes
Dpaixupipes=sortrows(CONDUITS,5);
%鈶℃柇i鏍规椂锛屾渶灏忓奖鍝嶇粍鍚堬細鍦―paixupipes閲岄?夊墠i鏍癸紱鏅?氶殢鏈虹粍鍚堬細鍦―paixupipes閲岄殢鏈洪?塱鏍癸紱鏈?澶у奖鍝嶇粍鍚堬細鍦―paixupipes閲岄?夊悗i鏍?
% generalzuhe=cell(1,200);
badRoughness=100;%鍧忔浖瀹?=100

lingdianerwu=Dpaixupipes(find(Dpaixupipes(:,5)==0.25),:);
for i=1:size(lingdianerwu,1)
    for j=1:size(uselesspipes,1)
        if lingdianerwu(i,1)==uselesspipes(j,1)
            wuyongpipes(i,:)=lingdianerwu(i,:);
        end
    end
end
wuyongpipes(find(wuyongpipes(:,1)==0),:)=[];
youyongpipes=setdiff(lingdianerwu,wuyongpipes,'rows');
xiaopaixu=[wuyongpipes;youyongpipes];
Dpaixupipes(1:size(lingdianerwu,1),:)=xiaopaixu;

%把出口管放到最后
Outpipes=dlmread('Outpipes.txt');
chukoupipes=zeros(size(Dpaixupipes,1),size(Dpaixupipes,2));
for i=1:size(Dpaixupipes,1)
    for j=1:size(Outpipes,2)
        if Dpaixupipes(i,1)==Outpipes(1,j)
            chukoupipes(i,:)=Dpaixupipes(i,:);
        end
    end
end
chukoupipes(find(chukoupipes(:,1)==0),:)=[];
for i=1:size(Dpaixupipes,1)
    for j=1:size(Outpipes,2)
        if Dpaixupipes(i,1)==Outpipes(1,j)
            Dpaixupipes(i,:)=0;
        end
    end
end
Dpaixupipes(find(Dpaixupipes(:,1)==0),:)=[];
Dpaixupipes=[Dpaixupipes;chukoupipes];


% ttt=size(Dpaixupipes,1);
LineToEdit=LineToEdit2;
Line=Line2;
Roughness=0.01;
% for r=1:200  %200绉嶆櫘閫氶殢鏈虹粍鍚堬紝200鍒楋紝姣忓垪鏄竴绉嶇粍鍚堬紙i鏍癸級
generalzuhe=zeros(1,size(Dpaixupipes,2));
% for e=1:size(Dpaixupipes,1) %e浠ｈ〃鎶藉嚑鏍?
dianxingbaoguanshu=[21 38 52];
Inputfile3=Inputfile;
lifetime=360;
for n=1:size(dianxingbaoguanshu,2) %e浠ｈ〃鎶藉嚑鏍?
    e=dianxingbaoguanshu(n);
%     if e>=2
%         while(1)
%             randompipe=Dpaixupipes(randperm(size(Dpaixupipes,1),1),:);%200绉嶆櫘閫氶殢鏈虹粍鍚堬紝200涓厓鑳烇紝姣忎釜鍏冭優鍖呭惈涓?绉嶇粍鍚堬紙i鏍癸級
%             repeat=find(generalzuhe(:,1)==randompipe(1,1));
%             if isempty(repeat)==1
%                 break
%             end    
%         end
%     else
%         randompipe=Dpaixupipes(randperm(size(Dpaixupipes,1),1),:);%200绉嶆櫘閫氶殢鏈虹粍鍚堬紝200涓厓鑳烇紝姣忎釜鍏冭優鍖呭惈涓?绉嶇粍鍚堬紙i鏍癸級
%     end
    randompipe=Dpaixupipes(end-(e-1):end,:);%200绉嶆櫘閫氶殢鏈虹粍鍚堬紝200涓厓鑳烇紝姣忎釜鍏冭優鍖呭惈涓?绉嶇粍鍚堬紙i鏍癸級
    generalzuhe(1:e,:)=randompipe;
    for i=1:size(CONDUITS,1) %姝ゅ惊鐜腑minzuhe鐨勯偅浜涚鐨凴oughness=badRoughness=100
        Temp=[CONDUITS(i,1)	CONDUITS(i,2)	CONDUITS(i,3)	CONDUITS(i,4)	Roughness	CONDUITS(i,9)	CONDUITS(i,8)	CONDUITS(i,10)	CONDUITS(i,11)];
        for k=1:size(generalzuhe,1)
            if CONDUITS(i,1)==generalzuhe(k,1)
                Temp=[CONDUITS(i,1)	CONDUITS(i,2)	CONDUITS(i,3)	CONDUITS(i,4)	badRoughness	CONDUITS(i,9)	CONDUITS(i,8)	CONDUITS(i,10)	CONDUITS(i,11)];
            end
        end
%         w(i,:)=Temp
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
    

    % [LID_USAGE]
    % ;;Subcatchment   LID Process      Number  Area       Width      InitSat    FromImp    ToPerv     RptFile                  DrainTo         
    % ;;-------------- ---------------- ------- ---------- ---------- ---------- ---------- ---------- ------------------------ ----------------
    %LID1 IT AREA=20 With=2 FromImp=45  LID2 RB AREA=1.5  FromImp=35
    %717    1	15	20	2	0	45	0 
    InitSat=0;
    AreaBC=20; %姣忎釜BC鍗曞厓200銕?
    WidthBC=10;
    AreaPP=15; %姣忎釜PP鍗曞厓150銕?
    WidthPP=10;
    ToPervBC=1; %鎵?鏈夋孩娴佽繑鍥炶嚦閫忔按鍖哄煙
    ToPervPP=1; %鎵?鏈夋孩娴佽繑鍥炶嚦閫忔按鍖哄煙
    LIDProcessBC=1;
    LIDProcessPP=2;

    myformat = '%2s%23s%6s%17s%11s%11s%11s%s%[^\n\r]';
    ToSearch=regexp( fileread('[LID_USAGE]-T.txt'), '\n', 'split');
    tf=strcmp(Inputfile,ToSearch);
    Line=find(tf);
    LineToEdit=Line+3;
    CapBC=zeros(1,size(Initial,1));
    CapPP=zeros(1,size(Initial,1));
    AreaofBC=zeros(1,size(Initial,1));
    AreaofPP=zeros(1,size(Initial,1));
    
    Inputfile2=Inputfile;
    LineToEdit3=LineToEdit;
    %PP-GREI
    for i=1:size(Initial,1)
    %      if dinp(i+119*2)==1 %1有LID
    %      if pardesign(i+size(pipes,1))>.5
         Subcatchment=Initial(i,1);
         usablearea = dlmread('usable area.txt');
         Area=usablearea(i,1);%汇水区面积%汇水区面积
    %          imperv=Initial(i,5);%不渗透比
         %Areaimperv=imperv*0.01*Area*10000;%不渗透部分面积
    %          Areaperv=(1-imperv*0.01)*Area*10000;%渗透部分面积
        if dinp(1,83*2+i*2)==0 && dinp(1,83*2+i*2-1)==0 %00 0%
           LIDSizeBC=0;
        elseif dinp(1,83*2+i*2)==1 && dinp(1,83*2+i*2-1)==0 %01 3.33%
           LIDSizeBC=Area*50/3;
        elseif dinp(1,83*2+i*2)==0 && dinp(1,83*2+i*2-1)==1 %10 6.67%
           LIDSizeBC=Area*50*2/3;
        elseif dinp(1,83*2+i*2)==1 && dinp(1,83*2+i*2-1)==1 %11 10%
           LIDSizeBC=Area*50;
        end
        if dinp(1,83*2+24+24+i*2)==0 && dinp(1,83*2+24+24+i*2-1)==0 %00 0%
           LIDSizePP=0;
        elseif dinp(1,83*2+24+24+i*2)==1 && dinp(1,83*2+24+24+i*2-1)==0 %01 3.33%
           LIDSizePP=Area*200/3/3;
        elseif dinp(1,83*2+24+24+i*2)==0 && dinp(1,83*2+24+24+i*2-1)==1 %10 6.67%
           LIDSizePP=Area*200/3*2/3;
        elseif dinp(1,83*2+24+24+i*2)==1 && dinp(1,83*2+24+24+i*2-1)==1 %11 10%
           LIDSizePP=Area*200/3;
        end
        if LIDSizeBC~=0 || LIDSizePP~=0
             nBC=round(LIDSizeBC); %nIT=round(Area*NumofITperHec) %子集水区面积*每公顷子集水区的BC比
             nPP=round(LIDSizePP); %nRB=round(Area*NumofRBperHec) %子集水区面积*每公顷子集水区的PP比
             %BC+PP面积不能大于该汇水区的10%
             BCbili=AreaBC*nBC/(AreaBC*nBC+AreaPP*nPP);
             PPbili=AreaPP*nPP/(AreaBC*nBC+AreaPP*nPP);
             if nBC*AreaBC+nPP*AreaPP>Area*10000 %如果加起来大于10%，按10%乘以BCPP比例
                 nBC=round((Area*10000*BCbili)/AreaBC);
                 nPP=round((Area*10000*PPbili)/AreaPP);
                 BCbili=AreaBC*nBC/(AreaBC*nBC+AreaPP*nPP);
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
%              Temp1=[Subcatchment	LIDProcessBC	nBC	AreaBC	WidthBC InitSat FromImpBC ToPervBC]; %Temp1=[Subcatchment	LIDProcessIT	nIT	AreaIT	WidthIT InitSat FromImpIT ToPervIT]; 
             Temp2=[Subcatchment	LIDProcessPP	nPP	AreaPP	WidthPP InitSat FromImpPP ToPervPP]; %Temp2=[Subcatchment	LIDProcessRB	nRB	AreaRB	WidthRB InitSat FromImpRB ToPervRB];
%              Temp1=num2str(Temp1);
             Temp2=num2str(Temp2);
%              Inputfile{1,LineToEdit} = sprintf(myformat,Temp1);
             Inputfile{1,LineToEdit} = sprintf(myformat,Temp2);
             LineToEdit=LineToEdit+1;
             AreaofBC(1,i)=nBC*AreaBC; %列出每个子汇水区的BC面积（公顷）
             CapBC(1,i)=75.6*AreaofBC(1,i)+15*(AreaofBC(1,i)^0.5); %126938.47
             AreaofPP(1,i)=nPP*AreaPP;
             CapPP(1,i)=34.6*AreaofPP(1,i)+15*(AreaofPP(1,i)^0.5); %126938.47
    %      else
         end
    end
    LineToEdit=51;
    Inputfile{1,LineToEdit} = sprintf(myformat,'1                VOLUME 0:01     1.0      TIMESERIES 5yChicago-6h');
    fid = fopen('LIMatF5PP.inp', 'w');
    fprintf(fid, '%s\n', Inputfile{:});
    fclose(fid);
    ! swmm5.exe  LIMatF5PP.inp  LIMatF5PP.rpt 
    filename = 'LIMatF5PP.rpt';
    delimiter = {'..','\t',',',' ',';'};
    startRow = 58;
    endRow = 58;
    formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    LIMatF = [dataArray{1:end-1}];
    a=LIMatF(1,1);
    if a=='Flooding'
    geneeralFloodvolume5yPP(e)=str2num(LIMatF(1,4));
    else 
    geneeralFloodvolume5yPP(e)=999;    
    end  
    VTFPP(e)=geneeralFloodvolume5yPP(e);

    %BC-GREI
    LineToEdit=LineToEdit3;
    Inputfile=Inputfile2;
    for i=1:size(Initial,1)
         Subcatchment=Initial(i,1);
         usablearea = dlmread('usable area.txt');
         Area=usablearea(i,1);%汇水区面积%汇水区面积
         if dinp(1,83*2+i*2)==0 && dinp(1,83*2+i*2-1)==0 %00 0%
           LIDSizeBC=0;
        elseif dinp(1,83*2+i*2)==1 && dinp(1,83*2+i*2-1)==0 %01 3.33%
           LIDSizeBC=Area*50/3;
        elseif dinp(1,83*2+i*2)==0 && dinp(1,83*2+i*2-1)==1 %10 6.67%
           LIDSizeBC=Area*50*2/3;
        elseif dinp(1,83*2+i*2)==1 && dinp(1,83*2+i*2-1)==1 %11 10%
           LIDSizeBC=Area*50;
        end
        if dinp(1,83*2+24+24+i*2)==0 && dinp(1,83*2+24+24+i*2-1)==0 %00 0%
           LIDSizePP=0;
        elseif dinp(1,83*2+24+24+i*2)==1 && dinp(1,83*2+24+24+i*2-1)==0 %01 3.33%
           LIDSizePP=Area*200/3/3;
        elseif dinp(1,83*2+24+24+i*2)==0 && dinp(1,83*2+24+24+i*2-1)==1 %10 6.67%
           LIDSizePP=Area*200/3*2/3;
        elseif dinp(1,83*2+24+24+i*2)==1 && dinp(1,83*2+24+24+i*2-1)==1 %11 10%
           LIDSizePP=Area*200/3;
        end
        if LIDSizeBC~=0 || LIDSizePP~=0
             nBC=round(LIDSizeBC); %nIT=round(Area*NumofITperHec) %子集水区面积*每公顷子集水区的BC比
             nPP=round(LIDSizePP); %nRB=round(Area*NumofRBperHec) %子集水区面积*每公顷子集水区的PP比
             %BC+PP面积不能大于该汇水区的10%
             BCbili=AreaBC*nBC/(AreaBC*nBC+AreaPP*nPP);
             PPbili=AreaPP*nPP/(AreaBC*nBC+AreaPP*nPP);
             if nBC*AreaBC+nPP*AreaPP>Area*10000 %如果加起来大于10%，按10%乘以BCPP比例
                 nBC=round((Area*10000*BCbili)/AreaBC);
                 nPP=round((Area*10000*PPbili)/AreaPP);
                 BCbili=AreaBC*nBC/(AreaBC*nBC+AreaPP*nPP);
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
%              Temp2=[Subcatchment	LIDProcessPP	nPP	AreaPP	WidthPP InitSat FromImpPP ToPervPP]; %Temp2=[Subcatchment	LIDProcessRB	nRB	AreaRB	WidthRB InitSat FromImpRB ToPervRB];
             Temp1=num2str(Temp1);
%              Temp2=num2str(Temp2);
             Inputfile{1,LineToEdit} = sprintf(myformat,Temp1);
%              Inputfile{1,LineToEdit} = sprintf(myformat,Temp2);
             LineToEdit=LineToEdit+1;
             AreaofBC(1,i)=nBC*AreaBC; %列出每个子汇水区的BC面积（公顷）
             CapBC(1,i)=75.6*AreaofBC(1,i)+15*(AreaofBC(1,i)^0.5); %126938.47
             AreaofPP(1,i)=nPP*AreaPP;
             CapPP(1,i)=34.6*AreaofPP(1,i)+15*(AreaofPP(1,i)^0.5); %126938.47
    %      else
         end
    end
    LineToEdit=51;
    Inputfile{1,LineToEdit} = sprintf(myformat,'1                VOLUME 0:01     1.0      TIMESERIES 5yChicago-6h');
    fid = fopen('LIMatF5BC.inp', 'w');
    fprintf(fid, '%s\n', Inputfile{:});
    fclose(fid);
    ! swmm5.exe  LIMatF5BC.inp  LIMatF5BC.rpt 
    filename = 'LIMatF5BC.rpt';
    delimiter = {'..','\t',',',' ',';'};
    startRow = 59;
    endRow = 59;
    formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    LIMatF = [dataArray{1:end-1}];
    a=LIMatF(1,1);
    if a=='Flooding'
    geneeralFloodvolume5yBC(e)=str2num(LIMatF(1,4));
    else 
    geneeralFloodvolume5yBC(e)=999;    
    end  
    VTFBC(e)=geneeralFloodvolume5yBC(e);
%         Tf(r)=floodduration5y(r);
%         Tn=16;
%     Resilience5yBC(e)=1-(VTFBC(e)/VPRBC(e)/1.0182);%total area=101.82ha %无LTE的LID-GREI
%         Resilience5y(r)=1-(VTF(r)/VTI(r))*(Tf(r)/Tn);

    %BCPP-GREI
    LineToEdit=LineToEdit3;
    Inputfile=Inputfile2;
    for i=1:size(Initial,1)
    %      if dinp(i+119*2)==1 %1有LID
    %      if pardesign(i+size(pipes,1))>.5
         Subcatchment=Initial(i,1);
         usablearea = dlmread('usable area.txt');
         Area=usablearea(i,1);%汇水区面积%汇水区面积
    %          imperv=Initial(i,5);%不渗透比
         %Areaimperv=imperv*0.01*Area*10000;%不渗透部分面积
    %          Areaperv=(1-imperv*0.01)*Area*10000;%渗透部分面积
        if dinp(1,83*2+i*2)==0 && dinp(1,83*2+i*2-1)==0 %00 0%
           LIDSizeBC=0;
        elseif dinp(1,83*2+i*2)==1 && dinp(1,83*2+i*2-1)==0 %01 3.33%
           LIDSizeBC=Area*50/3;
        elseif dinp(1,83*2+i*2)==0 && dinp(1,83*2+i*2-1)==1 %10 6.67%
           LIDSizeBC=Area*50*2/3;
        elseif dinp(1,83*2+i*2)==1 && dinp(1,83*2+i*2-1)==1 %11 10%
           LIDSizeBC=Area*50;
        end
        if dinp(1,83*2+24+24+i*2)==0 && dinp(1,83*2+24+24+i*2-1)==0 %00 0%
           LIDSizePP=0;
        elseif dinp(1,83*2+24+24+i*2)==1 && dinp(1,83*2+24+24+i*2-1)==0 %01 3.33%
           LIDSizePP=Area*200/3/3;
        elseif dinp(1,83*2+24+24+i*2)==0 && dinp(1,83*2+24+24+i*2-1)==1 %10 6.67%
           LIDSizePP=Area*200/3*2/3;
        elseif dinp(1,83*2+24+24+i*2)==1 && dinp(1,83*2+24+24+i*2-1)==1 %11 10%
           LIDSizePP=Area*200/3;
        end
        if LIDSizeBC~=0 || LIDSizePP~=0
             nBC=round(LIDSizeBC); %nIT=round(Area*NumofITperHec) %子集水区面积*每公顷子集水区的BC比
             nPP=round(LIDSizePP); %nRB=round(Area*NumofRBperHec) %子集水区面积*每公顷子集水区的PP比
             %BC+PP面积不能大于该汇水区的10%
             BCbili=AreaBC*nBC/(AreaBC*nBC+AreaPP*nPP);
             PPbili=AreaPP*nPP/(AreaBC*nBC+AreaPP*nPP);
             if nBC*AreaBC+nPP*AreaPP>Area*10000 %如果加起来大于10%，按10%乘以BCPP比例
                 nBC=round((Area*10000*BCbili)/AreaBC);
                 nPP=round((Area*10000*PPbili)/AreaPP);
                 BCbili=AreaBC*nBC/(AreaBC*nBC+AreaPP*nPP);
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
             AreaofBC(1,i)=nBC*AreaBC; %列出每个子汇水区的BC面积（公顷）
             CapBC(1,i)=75.6*AreaofBC(1,i)+15*(AreaofBC(1,i)^0.5); %126938.47
             AreaofPP(1,i)=nPP*AreaPP;
             CapPP(1,i)=34.6*AreaofPP(1,i)+15*(AreaofPP(1,i)^0.5); %126938.47
    %      else
         end
    end
    %
    LineToEdit=51;
    Inputfile{1,LineToEdit} = sprintf(myformat,'1                VOLUME 0:01     1.0      TIMESERIES 5yChicago-6h');%5yChicago-6h
    fid = fopen('LIMatF5.inp', 'w');
    fprintf(fid, '%s\n', Inputfile{:});
    fclose(fid);

    %LineToEdit鍜孡ine鍥炲埌鏈?鍒濈殑璧风偣
    LineToEdit=LineToEdit2;
    Line=Line2;
    Inputfile=Inputfile3;

    % RUN SWMM and read data 5y
    ! swmm5.exe  LIMatF5.inp  LIMatF5.rpt 
    filename = 'LIMatF5.rpt';
    delimiter = {'..','\t',',',' ',';'};
    startRow = 59;
    endRow = 59;
    formatSpec = '%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    LIMatF = [dataArray{1:end-1}];
    a=LIMatF(1,1);
    if a=='Flooding'
    geneeralFloodvolume5y(e)=str2num(LIMatF(1,4));
    else 
    geneeralFloodvolume5y(e)=999;    
    end  

    %Resilience calculation
    filename = 'LIMatF5.rpt';
    delimiter = {'..','\t',',',' ',';'};
    startRow = 42;
    endRow = 42;
    formatSpec = '%s%s%s%s%s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
%         formatSpec = '%*s%*s%*s%*s%*s%f%f%*s%*s%*s%*s%*s%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    LIMatF = [dataArray{1:end-2}];
    VPR(e)=str2num(LIMatF(1,4));
    VTF(e)=geneeralFloodvolume5y(e);
%         Tf(r)=floodduration5y(r);
%         Tn=16;
%         Resilience5y(r)=1-(VTF(r)/VTI(r))*(Tf(r)/Tn);
   
    % 加入BC性能衰退曲线
    baseBCcontrol=VTFPP(e)-VTF(e);%PP+GREI的flooding - BC+PP+GREI的flooding = BC的基准控制量
    pdfnorm_BCstan=zeros(1,13);
    for x=1:13
        pdfnorm_BCstan(x)=(1/(2*pi)^0.5)*exp(1)^(-(x-7)^2/2)+1-((1/(2*pi)^0.5)*exp(1)^(-0^2/2));
    end
    %12个月的BC性能衰退曲线
    pdfnorm_BC=zeros(1,lifetime+1);
    pdfnorm_BC(1:12+1)=pdfnorm_BCstan;
    for i=2:lifetime/12
        pdfnorm_BC(12*(i-1)+2:12*i+1)=pdfnorm_BCstan(2:13);
    end
    % plot(1:361,pdfnorm_BC);%600个月的BC性能（不衰退）曲线
    pdfnorm_BC0=1;
    LSEhighest_BC=zeros(1,lifetime/12+1);
    for i=1:lifetime/12+1
        LSEhighest_BC(i)=(1-0.02)^(i-1);%BC性能最高点衰退曲线
    end
    LSEmean_BC=zeros(1,lifetime+1);
    for i=1:lifetime+1
        LSEmean_BC(i)=pdfnorm_BC(i)*LSEhighest_BC(fix((i-2)/12)+1)/pdfnorm_BC0;
    end

    % 拟合BC性能衰退曲线,注释则结果呈波动
    fit=polyfit(0:lifetime,LSEmean_BC,2); %拟合2项式
    x=(0:lifetime);
    LSEmean_BC=fit(1)*x.^2+fit(2)*x+fit(3);
    BCcontrolofmonth=baseBCcontrol*LSEmean_BC(120*n+1);%考虑衰退后的BC控制量
    
    % 加入PP性能衰退曲线
    basePPcontrol=VTFBC(e)-VTF(e);%BC+GREI的flooding - BC+PP+GREI的flooding = PP的基准控制量
    LSEmean_PP=zeros(1,lifetime+1);
    m=1:(1/12):lifetime/12+1;
    for i=1:lifetime+1
        LSEmean_PP(i)=-0.02*(m(i)-1)+1;%PP性能衰退曲线
    end
    %PP性能衰退曲线
    PPcontrolofmonth=basePPcontrol*LSEmean_PP(120*n+1);%考虑衰退后的PP控制量
    FVbaseLID_GERIwithoutLID=VTF(e)+basePPcontrol+baseBCcontrol;%无LID的LID-GREI
    %%%

    LID_GREIofmonth=FVbaseLID_GERIwithoutLID-PPcontrolofmonth-BCcontrolofmonth;%考虑衰退后的LID-GREI=无LID的LID-GREI-衰退后的PP控制量-衰退后的BC控制量
    Resilience5y(e)=1-(LID_GREIofmonth/VPR(e)/1.0182);%total area=101.82ha %无LTE的LID-GREI

end        
%e=2-118鏃舵瘡涓猠鐨刦inalFloodvolume
% geneeralFloodvolume5y(geneeralFloodvolume5y==0)=[];
Resilience5y(Resilience5y==0)=[];
% maxFloodvolume5y=geneeralFloodvolume5y;
maxResilience5y=Resilience5y;
% end
Res_OPR=(maxResilience5y(1)*maxResilience5y(2)*maxResilience5y(3))^(1/3);
