clc
clear
global Element GaussXW Node LE MatiA Mat AG  ElNoC DoF
[ProblemType, DoF,Mat,Node,Area,Element,BCU,InNode]=Input('data110.txt');
AG=[0 -9.8]';       %   acceleration of gravity
IntMethod=3;      %0:解析解   1:int    2:quadv   3:高斯积分
GaussPt=20;
[x, w] = GaussLegendre(GaussPt);
GaussXW = [x w]';
 if ProblemType==2
     for  ia=1:size(Mat,1)
       %  Mat(ia,2)=Mat(ia,2)*(1+2*Mat(ia,3))/(1+Mat(ia,3))^2;
         Mat(ia,3)=Mat(ia,3)/(1+Mat(ia,3));
    end 
 end
MChau=0;
% 以下是程序部分
set(0,'FormatSpacing','compact')
%前处理中公共边界的处理
%公共边界的查找
CommonElement=zeros(0,4);
for iE=1:size(Element,1)
    IP=find(Element((iE+1):end,4)==Element(iE,3) &...
        Element((iE+1):end,3)==Element(iE,4));
    if ~isempty(IP)
       iE1=Element(iE,2);
       MatA1=Element(iE,1);
       MatA2=Element(iE+IP,1);
       if MatA1<MatA2
           CommonElement=cat(1,CommonElement,[iE1 MatA1 0 MatA2]);
       else
           CommonElement=cat(1,CommonElement,[iE1 MatA2 0 MatA1]);
       end
    end
end
%增加节点，修改对应单元中的节点号并增加其对应的边界条件

Element(:,end+2)=0;
for iA0=1:(size(Area,1)-1)
    for iA1=(iA0+1):size(Area,1)
        IndexComEl=find(CommonElement(:,2)==iA0 & CommonElement(:,4)==iA1);
        if isempty(IndexComEl)
            continue
        end
        for iCE=1:length(IndexComEl)
            E0=CommonElement(IndexComEl(iCE),1);
            E1=max(Element(:,2))+1;
            InE1=find(Element(:,2)==CommonElement(IndexComEl(iCE),1) & ...
            Element(:,1)==CommonElement(IndexComEl(iCE),4));
            Element(InE1,2)=E1;
            CommonElement(IndexComEl(iCE),3)=E1;
            InE0=find(Element(:,2)==E0);
            for eN = 1:2
                N0=Element(InE0,eN+2);
                N1=Element(InE1,5-eN);
%                李春光的程序
                inENa=find(Element(:,3)==N1 & Element(:,1)>iA0 & Element(:,13)==0);
                inENb=find(Element(:,4)==N1 & Element(:,1)>iA0 & Element(:,14)==0);
%                 邓琴的程序
%                 temp=iA0;
%                 if Element(InE0,13:14)==0 & Element(InE1,13:14)~=0
%                     temp=Element(InE1,1);
%                 end
%                 inENa=find(Element(:,3)==N1 & Element(:,1)~=temp & Element(:,13)==0);
%                 inENb=find(Element(:,4)==N1 & Element(:,1)~=temp & Element(:,14)==0);
                if  isempty(union(inENa,inENb)) & N0~=N1
                    continue
                end
                N11=size(Node,1)+1;
                
                if isempty(union(inENa,inENb)) %前面大面号，后面小面号
                    Node(N11,:)=[N11 Node(N0,2:3)];
                    inENa=find(Element(:,3)==N0 & Element(:,1)<iA1 & Element(:,13)==0);
                    inENb=find(Element(:,4)==N0 & Element(:,1)<iA1 & Element(:,14)==0);
                    Element(inENa,3)=N11;
                    Element(inENb,4)=N11;
                    Element(find(Element(:,3)==N11 & Element(:,1)==iA0),13)=1;
                    Element(find(Element(:,4)==N11 & Element(:,1)==iA0),14)=1;
                else    %前面小面号，后面大面号
                    Node(N11,:)=[N11 Node(N1,2:3)];
                    Element(inENa,3)=N11;
                    Element(inENb,4)=N11;
                    Element(find(Element(:,3)==N11 & Element(:,1)==iA1),13)=1;
                    Element(find(Element(:,4)==N11 & Element(:,1)==iA1),14)=1;
                end                
            end
            Element(InE0,5:8)=0;
            Element(InE1,5:8)=0;
        end        
    end
end
Element(:,13:end)=[];
%对单元进行排序
 [B,IX]=sort(Element(:,2));
 Element=Element(IX,:);
 
CommonNodes=zeros(0,2);  %[主节点 从节点]
for iE=1:size(CommonElement,1)
    for eN=1:2
        Na=Element(CommonElement(iE,1),eN+2);
        Nb=Element(CommonElement(iE,3),5-eN);
        if isempty(find(CommonNodes(:,1)==Na & CommonNodes(:,2)==Nb))
            if isempty(find(CommonNodes(:,1)==Na))  %每个节点在主节点只出现一次
                CommonNodes=cat(1,CommonNodes,[Na Nb]);
            else
                CommonNodes=cat(1,CommonNodes,[Nb Na]);
            end
        end
    end
end
%去除多余的公共点对,指该点四周全被各个面包围的点，不去除会增加多余的约束
RowTobeDel=zeros(0,1);
CommonNodes(:,end+1)=0;
for iE=1:size(CommonNodes,1)
    UsedRow=[];
    if CommonNodes(iE,end)==1
        continue
    end
    Node0=CommonNodes(iE,1);
    Nodei=CommonNodes(iE,2);
    Nodej=0;
    Row = find(CommonNodes(:,1)==Nodei);
    if CommonNodes(Row,2)~=Node0
        Nodej=CommonNodes(Row,2);
        UsedRow=cat(1,UsedRow,Row);
    end
    Row = find(CommonNodes(:,2)==Nodei);
    if CommonNodes(Row,1)~=Node0
        Nodej=CommonNodes(Row,1);
        UsedRow=cat(1,UsedRow,Row);
    end
    while  Nodej~=0        
        Row=find((CommonNodes(:,1)==Nodej) & (CommonNodes(:,2)~=Nodei)...
                |(CommonNodes(:,2)==Nodej) & (CommonNodes(:,1)~=Nodei));
        if isempty(Row)
            break
        end
        CommonNodes(Row,end)=1;
        Nodei=Nodej;
        Nodej=sum(CommonNodes(Row,1:2))-Nodei;
        if Nodej==Node0
            RowTobeDel=cat(1,RowTobeDel,Row);
            break
        end            
    end 
    CommonNodes(iE,end)=1;
    CommonNodes(UsedRow,end)=1;
end
CommonNodes(RowTobeDel,:)=[];
CommonNodes(:,end)=[];
% 求单元的法向
 DL=Node(Element(:,4),2:3)-Node(Element(:,3),2:3);
 LE=sum(DL.*DL,2).^0.5;
 DL=DL./[LE LE];
 NE=DL*[0 -1;1 0];
%边界节点已知位移量个数
 UPknow=zeros(size(Node,1),7);
 UPknow(:,1)=1:size(Node,1);
 IndexUPknow=zeros(size(Node,1),8);
 IndexUPknow(:,1)=1:size(Node,1);
 IUX=find(BCU(:,2)==1);
 IUY=find(BCU(:,2)==2);
%  IUXY=find(BCU(:,2)==3);
 
 UPknow(BCU(IUX),2)=BCU(IUX,3);
 UPknow(BCU(IUY),3)=BCU(IUY,3);
%  UPknow(BCU(IUXY),2:3)=BCU(IUXY,3)*[1 1]

 IndexUPknow(BCU(IUX),2)=1;
 IndexUPknow(BCU(IUY),3)=1;
%  IndexUPknow(BCU(IUXY),2:3)=ones(size(IUXY,1),2)
 
 
IPXb=Element(find(Element(:,7)==1),4);
IPYb=Element(find(Element(:,8)==1),4);
% IPXYb=find(BCP(:,2)==3)
IPXa=Element(find(Element(:,5)==1),3);
IPYa=Element(find(Element(:,6)==1),3);
% IPXYa=find(BCP(:,2)==13)

UPknow(IPXb,4)=Element(find(Element(:,7)==1),11);
UPknow(IPYb,5)=Element(find(Element(:,8)==1),12);
% UPknow(BCP(IPXYb),4:5)=BCP(IPXYb,3)*[1 1]
UPknow(IPXa,6)=Element(find(Element(:,5)==1),9);
UPknow(IPYa,7)=Element(find(Element(:,6)==1),10);
% UPknow(BCP(IPXYa),6:7)=BCP(IPXYa,3)*[1 1]

IndexUPknow(IPXb,4)=1;
IndexUPknow(IPYb,5)=1;
% IndexUPknow(BCP(IPXYb),4:5)=ones(size(BCP(IPXYb,3)))*[1 1]
IndexUPknow(IPXa,6)=1;
IndexUPknow(IPYa,7)=1;
% IndexUPknow(BCP(IPXYa),6:7)=ones(size(BCP(IPXYa,3)))*[1 1]

IndexUPknow(:,end)=sum(IndexUPknow(:,2:7),2);
IndexUPknow(:,end+1)=IndexUPknow(:,end);

%IndexU1已知位移边界，IndexU0未知位移边界

U1=cat(1,[BCU(IUX,1)*2-1,BCU(IUX,3)],[BCU(IUY,1)*2,BCU(IUY,3)]);
[IndexU1,IX] = sort(U1(:,1));
U1=U1(IX,2);
IndexU0=setdiff(1:size(Node,1)*2, IndexU1)';

%IndexP1已知应力边界，IndexP0未知应力边界
 P1=cat(1,[IPXb*4-3,UPknow(IPXb,4)],[IPYb*4-2,UPknow(IPYb,5)],...
[IPXa*4-1,UPknow(IPXa,6)],[IPYa*4,UPknow(IPYa,7)]);
 [IndexP1,IX] = sort(P1(:,1));
 P1=P1(IX,2);
IndexP0=setdiff(1:size(Node,1)*4, IndexP1)';

tic


for iA=1:max(Element(:,1))
    MatiA=Area(iA,2);
    Den=Mat(MatiA,1);
    NU=Mat(MatiA,3);
    SG=Mat(MatiA,4);
    ElementAi=Element(find(Element(:,1)==iA),:);
    if DoF==0
       [H,G,VF]=HGFA(ElementAi);
    end
end

Gm=mean(Mat(:,4));
G=G*Gm;
P1=P1/Gm;
%计算
HGL=cat(2,H(:,IndexU0),-G(:,IndexP0));
HGR=cat(2,-H(:,IndexU1),G(:,IndexP1));
UP1=cat(1,U1,P1);
F=HGR*UP1;
F=F+VF;

%施加位移连续条件
for eN=1:size(CommonNodes,1)
    iN0=CommonNodes(eN,1);    %主节点号
    iN1=CommonNodes(eN,2);    %从节点号
    for iDOF=1:2             
        ColPos1=numel(IndexUPknow(1:iN1-1,2:3))-sum(sum(IndexUPknow(1:iN1-1,2:3)))+...
            1+(1-IndexUPknow(iN1,2))*(iDOF-1);           
        if IndexUPknow(iN0,iDOF+1)==0        %位移未知
            ColPos0=numel(IndexUPknow(1:iN0-1,2:3))-sum(sum(IndexUPknow(1:iN0-1,2:3)))+...
            	1+(1-IndexUPknow(iN0,2))*(iDOF-1);
        	HGL(end+1,[ColPos0 ColPos1])=[1 -1];
            F(end+1)=0;
        else
            HGL(end+1,ColPos1)=1;
            F(end+1)=BCU(find(BCU(:,1)==iN0 & BCU(:,2)==iDOF),3);
        end
    end
    IndexUPknow(iN1,end)=IndexUPknow(iN1,end)+2; %所有位移连续条件给从节点
end
%施加剪应力正应力相等条件iDOF代表前后
for iCE=1:size(CommonElement,1)
    for eN=1:2
        iNM=Element(CommonElement(iCE,1),2+eN);
        iNS=Element(CommonElement(iCE,3),5-eN);
        ColPosM=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
            4*(iNM-1)-sum(sum(IndexUPknow(1:iNM-1,4:7)))+...
            (2-sum(IndexUPknow(iNM,4:5)))*(2-eN);
        ColPosS=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
            4*(iNS-1)-sum(sum(IndexUPknow(1:iNS-1,4:7)))+...
            (2-sum(IndexUPknow(iNS,4:5)))*(eN-1);  
        HGL(end+(1:2),[ColPosM+(1:2) ColPosS+(1:2)])=[eye(2) eye(2)];
        F(end+(1:2))=0;
        
        if CommonElement(iCE,1)==Element(find(Element(:,3)==iNM),2)
            E2M=Element(find(Element(:,4)==iNM),2);
        else
            E2M=Element(find(Element(:,3)==iNM),2);
        end
        E2S=CommonElement(find(CommonElement(:,1)==E2M),3);
        bMidNode=0;
        if ~isempty(E2S) & Element(E2S,1)==Element(CommonElement(iCE,3),1)
            bMidNode=1;%共边的中间节点
        end       
        
        if bMidNode==1 
%             |...
%            (IndexUPknow(iNM,end) < 4 & IndexUPknow(iNS,end) < 4)
            IndexUPknow(iNM,end)=IndexUPknow(iNM,end)+1;
            IndexUPknow(iNS,end)=IndexUPknow(iNS,end)+1;
        else
            for iDOF=1:2        %所有的连续条件首先满足方程少的点
%                 以前的程序
%                 if IndexUPknow(iNM,end)<= IndexUPknow(iNS,end)
%                     IndexUPknow(iNM,end)=IndexUPknow(iNM,end)+1; 
%                 else
%                     IndexUPknow(iNS,end)=IndexUPknow(iNS,end)+1;
%                 end
                Rows=[iNM iNS];
                [CMin IMin] = min(IndexUPknow(Rows,end));
                RMin=Rows(IMin);
                IndexUPknow(RMin,end)=IndexUPknow(RMin,end)+1;
            end            
        end
    end 
end
%分析、判断----补充方程
NodeNo=find(IndexUPknow(:,9)==3 | IndexUPknow(:,9)==2);
for eN=1:numel(NodeNo)
    Nodei=NodeNo(eN);    
    Eb=find(Element(:,3)==Nodei);
    Ea=find(Element(:,4)==Nodei);
    N1=NE(Eb,:);
    N2=NE(Ea,:);
    N21=[N2 -N1];
    PositonP=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
        4*(Nodei-1)-sum(sum(IndexUPknow(1:Nodei-1,4:7)));
    indexP0 = find(IndexUPknow(Nodei,4:7)==0);
    indexP1 = setdiff([1:4],indexP0);
    HGL(end+1,PositonP+(1:numel(indexP0)))=N21(indexP0);
    F(end+1)=-N21(indexP1)*UPknow(Nodei,3+indexP1)';

    if IndexUPknow(Nodei,end)==3 & norm(N1-N2)<1e-3
        if  IndexUPknow(Nodei,2)==0 & IndexUPknow(Nodei,3)==0
            HGL(end,:)=[];
            F(end)=[];
            continue
        else 
            continue
        end
     
    end
    IndexRow=size(HGL,1)+1;
    if  norm(N1-N2)<1e-3 %& IndexUPknow(Nodei,end)~=3 %两单元共线，特殊的角点约束
        MN=(N1+N2)/2;
        N21=[MN(2) -MN(1) -MN(2) MN(1)];
        HGL(IndexRow,PositonP+(1:numel(indexP0)))=N21(indexP0);
        F(IndexRow)=-N21(indexP1)*UPknow(Nodei,3+indexP1)';
    elseif MChau==1%加另一个角点约束            
        G2=2*Mat(Area(Element(Eb,1),2),4)/Gm;
        CChauN1=-G2/LE(Eb)*[-N1(2) N1(1)];
        CChauN2=-G2/LE(Ea)*[-N2(2) N2(1)];
        CChauN0=-CChauN1-CChauN2;
        CChauU=[CChauN1;CChauN0;CChauN2];
        CChauP=[N1 -N2];
        NChau=[Element(Eb,3) Nodei Element(Ea,4)];
        F(IndexRow)=0;
        for iNChau=1:3
            IndexColumn=2*(NChau(iNChau)-1)-sum(sum(IndexUPknow(1:(NChau(iNChau)-1),2:3)));
            for iDof=1:2
                if IndexUPknow(NChau(iNChau),iDof+1)==0 %自由度未知   %位移未知
                    IndexColumn=IndexColumn+1;
                    HGL(IndexRow,IndexColumn)=CChauU(iNChau,iDof);
                else    %自由度已知     %位移已知
                    UiDof=BCU(find(BCU(:,1)==NChau(iNChau) & BCU(:,2)==iDof),3);
                    F(IndexRow)=F(IndexRow)-CChauU(iNChau,iDof)*UiDof;
                end
            end
        end
        HGL(IndexRow,PositonP+(1:4-sum(IndexUPknow(Nodei,4:7))))=CChauP(indexP0);
        F(IndexRow)=F(IndexRow)-CChauP(indexP1)*UPknow(Nodei,3+indexP1)';
    else
        %-------------------虚节点法-------------------begin
        
        HGL(end,:)=[];
        F(end)=[];
        
        ENo=find(Element(:,3)==Nodei);
        Nodej=Element(ENo,4);
        Nodek=size(Node,1)+1;
        Node(end+1,:)=[size(Node,1)+1 sum(Node([Nodei Nodej] ,2:3))/2];
        
        Element(end+1,:)=Element(ENo,:);
        Element(ENo,4)=Nodek;
        Element(end,[2 3])=[size(Element,1) Nodek];
        LE(ENo)=LE(ENo)/2;
        LE(end+1)=LE(ENo);
        NE(end+1,:)=NE(ENo,:);        
        
        np=Nodek;
        iA=Element(ENo,1);
        ElementAi=Element(find(Element(:,1)==iA),:);
        NodeAi=ElementAi(:,3:4);
        disp(sprintf('Area=%d        iN=%d',iA,np))
        
        %求C矩阵 
        n1=ElementAi(find(ElementAi(:,4)==np),3);
        n2=ElementAi(find(ElementAi(:,3)==np),4);
        a1=DirectionAngle(Node(np,2:end),Node(n1,2:end));
        a2=DirectionAngle(Node(np,2:end),Node(n2,2:end));
        if a2<a1
            a2=a2+2*pi;
        end
        da=a2-a1;
        aa=(a1+a2);
        nu2=2*(1-NU);
        Ci=eye(2)-1/2/pi/nu2*[nu2*da+cos(aa)*sin(da) sin(aa)*sin(da);...
                              sin(aa)*sin(da) nu2*da-cos(aa)*sin(da)];
        RIndex=2*np+[-1 0];
        C(RIndex,RIndex)=Ci;
        
        xp=Node(np,2:end)';
        [Hnp Gnp VFnp]=HGFA(np,ElementAi,IntMethod);
        
        Hnp(:,end-1:end)=Hnp(:,end-1:end)+Ci;
        HalfOneN=0.5*[diag(ones(2,1)) diag(ones(2,1))];
        Hnp(:,[Nodei*2+(-1:0) Nodej*2+(-1:0)])=...
            Hnp(:,[Nodei*2+(-1:0) Nodej*2+(-1:0) Nodek*2+(-1:0)])*...
            [diag(ones(4,1));HalfOneN];
        Hnp(:,Nodek*2+(-1:0))=[];

        Gnp=Gnp*Gm;
        Gnp(:,[Nodei*4+(-1:0) Nodej*4+(-3:-2)])=...
            Gnp(:,[Nodei*4+(-1:0) Nodej*4+(-3:-2) Nodek*4+(-3:0)])*...
            [diag(ones(4,1));HalfOneN;HalfOneN];
        Gnp(:,Nodek*4+(-3:0))=[];
        
        %计算
        HGL(end+1:end+2,:)=cat(2,Hnp(:,IndexU0),-Gnp(:,IndexP0));
        HGRnp=cat(2,-Hnp(:,IndexU1),Gnp(:,IndexP1));
        HGR(end+1:end+2,:)=HGRnp;        
        F(end+1:end+2)=HGRnp*UP1+VFnp;
        

        
        Node(end,:)=[];
        Element(end,:)=[];
        Element(ENo,4)=Nodej;
        LE(ENo)=LE(ENo)*2;
        LE(end)=[];
        NE(end,:)=[];
        %-------------------虚节点法-------------------end
        
        %-------------------临单元剪应力约束法-------------------begin
%         Nodej=Element(find(Element(:,3)==Nodei),4);
%         Nodek=Element(find(Element(:,3)==Nodej),4);
%         PositonP=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
%             4*(Nodei-1)-sum(sum(IndexUPknow(1:Nodei-1,4:7)));
%         Posi=PositonP;
%         PositonP=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
%             4*(Nodej-1)-sum(sum(IndexUPknow(1:Nodej-1,4:7)));
%         Posj=PositonP;
%         PositonP=2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))+...
%             4*(Nodek-1)-sum(sum(IndexUPknow(1:Nodek-1,4:7)));
%         Posk=PositonP;
%         HGL(end+1,[Posi Posj Posk])=[1 -2 1];
%          F(end+1)=0;
%         HGL(end+1,[Posi Posj Posk]+1)=[1 -2 1];
%          F(end+1)=0;
         %-------------------临单元剪应力约束法-------------------end
    end 
end

%C------------------------------
UP0=HGL\F;
U0=UP0(1:(2*size(Node,1)-sum(sum(IndexUPknow(:,2:3)))));
P0=UP0((2*size(Node,1)-sum(sum(IndexUPknow(:,2:3))))+1:end);
P0=P0*Gm;
Utemp=zeros(size(Node,1),2)';
Utemp(find(IndexUPknow(:,2:3)'==0))=U0;
UPknow(:,2:3)=UPknow(:,2:3)+Utemp';
Ptemp=zeros(size(Node,1),4)';
Ptemp(find(IndexUPknow(:,4:7)'==0))=P0;
UPknow(:,4:7)=UPknow(:,4:7)+Ptemp';
disp('boundary point:');
disp(sprintf('%s\t\t\t\t%s\t\t\t\t\t\t\t%s','Node No','UX','UY'));
disp(sprintf('%d\t%24.12f\t%24.12f\n',UPknow(:,1:3)'));
disp(sprintf('%s\t\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s','Node No','SXb','SYb','SXa','SYa'))
disp(sprintf('%d\t%16.4f\t%16.4f\t%16.4f\t%16.4f\n',UPknow(:,[1 4:7])'))


% 计算内部点的位移
InCount=size(InNode,1);
iNterU=zeros(0,5);
iNterP=zeros(0,6);
for iKP=1:InCount
    xp=(InNode(iKP,1:2))';
    iKPState=-1;           %IN(any area)=1  ON=0  OHTERS=-1
    for iA=1:max(Element(:,1))
        MatiA=Area(iA,2);
        Den=Mat(MatiA,1);
        NU=Mat(MatiA,3);
        SG=Mat(MatiA,4);   
        EA=Element(find(Element(:,1)==iA),3:4);
        LoopNo=1;
        ON=0;
        while size(EA,1)~=0
        
            LoopNode=zeros(0,1);
            [Cmin Imin]=min(Node(EA(:,1),2));
            LoopNode(1:2,1)=EA(Imin,:)';
            EA(Imin,:)=[];
            while LoopNode(1)~=LoopNode(end)
                IRow=find(EA(:,1)==LoopNode(end));
                LoopNode(end+1)=EA(IRow,2);
                EA(IRow,:)=[];
            end
            LoopNode(end)=[];
            [INLoop ONLoop] = inpolygon(xp(1),xp(2),Node(LoopNode,2),Node(LoopNode,3));
            if ONLoop==1
                ON=1;
                break
            end            
            if LoopNo==1 
                if INLoop==0
                    IN=0;
                    break
                else
                    IN=1;
                end
            else
                if INLoop==1
                    IN=0;
                    break
                end
            end
            LoopNo=LoopNo+1;
        end
        
        ElementAi=Element(find(Element(:,1)==iA),:);        
         if ON==1
            ElementAiNode=[ElementAi(:,2) Node(ElementAi(:,3),2:3) Node(ElementAi(:,4),2:3)];
            iE=ElementAi(find(((ElementAiNode(:,2)-xp(1)).*(ElementAiNode(:,4)-xp(1))<=0).*...
            ((ElementAiNode(:,3)-xp(2)).*(ElementAiNode(:,5)-xp(2))<=0)),2);
            iE=iE(1);
            Nodeb=Element(iE,3);
            Nodea=Element(iE,4); 
            Lxpa=norm(xp'-Node(Nodea,2:3));
            Rb=Lxpa/LE(Element(iE,2));
            Ra=1-Rb;
            iNterUi=Ra*UPknow(Nodea,2:3)+Rb*UPknow(Nodeb,2:3);
            iNterU=cat(1,iNterU,[iA InNode(iKP,:) iNterUi]);
            
            PA=Ra*UPknow(Nodea,4:5)+Rb*UPknow(Nodeb,6:7);
            S11bar=PA*NE(iE,:)';
            S12bar=PA*DL(iE,:)';
            Strain22bar=diff(UPknow([Nodeb Nodea],2:3))*DL(iE,:)'/LE(iE);
            S22bar=1/(1-NU)*(NU*S11bar+2*SG*Strain22bar);
            
            MC=[NE(iE,1) NE(iE,2);-NE(iE,2) NE(iE,1)];
            Sxp=MC'*[S11bar S12bar;S12bar S22bar]*MC;

            iNterP=cat(1,iNterP,[iA InNode(iKP,:) [Sxp(1) Sxp(2) Sxp(4)]]);
            continue
       end  
       if IN==1
            AUi=zeros(2,1);
            APi=zeros(2,1);
            ADi=zeros(3,1);
            ASi=zeros(3,1);
            FVjj=zeros(2,1);
            SVjj=zeros(3,1);

        %   单元循环
            for iE=1:size(ElementAi,1)
                AD=zeros(3,4);
                AS=zeros(3,4);
                G=zeros(2,4);
                H=zeros(2,4);
                SV=zeros(3,1);
                SV11=0;
                SV12=0;
                SV22=0;
                FV=0;
                ElNo=ElementAi(iE,2);
                En1=Element(ElNo,3);
                En2=Element(ElNo,4);
                dpr1=xp-Node(En1,2:3)';dp1=sqrt(dpr1.'*dpr1);
                dpr2=xp-Node(En2,2:3)';dp2=sqrt(dpr2.'*dpr2);
                if dp1<LE(ElNo) && dp2<LE(ElNo)
                        [Gi,Hi,ADii,ASii,FVi,SVi]=OUTP(Node([En1 En2],:),Mat(MatiA,:),xp,AG,200);
                        G=G+Gi;
                        H=H+Hi;
                        AD=AD+ADii;
                        AS=AS+ASii;
                        FV=FV+FVi;
                        SV=SV+SVi;

                else
                      [Gi, Hi, ADii,ASii, FVi, SVi]=OUTP(Node([En1 En2],:),Mat(MatiA,:),xp,AG,GaussPt);
                      G=G+Gi;
                      H=H+Hi;
                      AD=AD+ADii;
                      AS=AS+ASii;
                      FV=FV+FVi;
                      SV=SV+SVi;
                end
                AUi(:,1)=AUi(:,1)+G*[UPknow(En1,6:7) UPknow(En2,4:5)]';
                APi(:,1)=APi(:,1)+H*[UPknow(En1,2:3) UPknow(En2,2:3)]';
                ADi(:,1)= ADi(:,1)+AD*[UPknow(En1,6:7) UPknow(En2,4:5)]';
                ASi(:,1)= ASi(:,1)+AS*[UPknow(En1,2:3) UPknow(En2,2:3)]';
                FVjj(:,1)=FVjj(:,1)+FV;
                SVjj(:,1)=SVjj(:,1)+SV;    
            end
            iNterUi=(AUi-APi+FVjj)';
            iNterPi=(ADi-ASi+SVjj)';
           iNterU=cat(1,iNterU,[iA InNode(iKP,:) iNterUi]);
           iNterP=cat(1,iNterP,[iA InNode(iKP,:) iNterPi]);
           break
       end

    end
end
if InCount~=0
    disp('internal point:');
    disp(sprintf('\t%s\t\t%s\t\t\t\t%s\t\t\t\t\t\t%s\t\t\t\t\t\t\t%s','Area','CX','CY','UX','UY'));
    disp(sprintf('%5.0f\t%12.5f\t%12.5f\t%24.12f\t%24.12f\n',iNterU'));
    disp(sprintf('\t%s\t\t%s\t\t\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s\t\t\t\t\t%s','Area','CX','CY','SX','SXY','SY'))
    disp(sprintf('%5.0f\t%12.5f\t%12.5f\t%16.4f\t%16.4f\t%16.4f\n',iNterP'))
end

toc
% close all
figure(1)
PlotSketch(0,UPknow(:,1:3))
PlotSketch(1,UPknow(:,1:3))