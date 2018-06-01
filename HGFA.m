function [Hnp,Gnp,VFnp]=HGFA(ElementAi)
global Element  Node LE MatiA Mat AG ElNoC DoF GaussXW
RowNos=size(Element,1)*2;
ColuNos=size(Element,1)*2;
if  DoF~=0
     RowNos=RowNos+DoF*2;
end
Hnp=zeros(RowNos);
Gnp=zeros(RowNos,ColuNos*2);
VFnp=zeros(ColuNos,1);
Rindex=1:2;
% VFnp=zeros(2,1);
for iA=1:size(ElementAi,1)
    disp(sprintf('Area=%d        iA=%d',iA,ElementAi(iA,2)))
    ElNoy=ElementAi(iA,2);   %外层积分单元
    En1y=Element(ElNoy,3);  %外层积分单元节点
    En2y=Element(ElNoy,4);

    for iW=1:size(GaussXW,2)
        x=GaussXW(1,iW);
        xp=[1-x x+1]/2*[Node(En1y,2:3);Node(En2y,2:3)];
        Hx=zeros(2,size(Node,1)*2);
        Px=zeros(2,4);
        Gx=zeros(2,size(Node,1)*4);
        VFx=zeros(2,1);
        for iE=1:size(ElementAi,1)
             ElNox=ElementAi(iE,2);
             En1x=Element(ElNox,3);
             En2x=Element(ElNox,4);
             if ElNox ~= ElNoy
                [Hi,Pi,Gi,VFi]=HGFE(xp',GaussXW,Node([En1x En2x],:),Mat(MatiA,:),AG);  %源点位于单元外部
             elseif ElNox == ElNoy
                [Hi,Pi,Gi,VFi]=HGCE(xp',Node([En1x En2x],:),Mat);   %源点位于单元内部
             end 
             CIndex=[2*En1x+[-1 0],2*En2x+[-1 0]];
             Hx(Rindex,CIndex)=Hx(Rindex,CIndex)+Hi;
             Px=Px+Pi*[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2];
             CIndex=[4*En1x+[-1 0] 4*En2x+[-3 -2]];
             Gx(Rindex,CIndex)=Gx(Rindex,CIndex)+Gi;
             VFx(Rindex)=VFx(Rindex)+VFi;
        end
        RIndex=[2*En1y+[-1 0],2*En2y+[-1 0]];
        Hnp(RIndex,:)=Hnp(RIndex,:)+[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2]'*Hx*GaussXW(2,iW)*LE(ElNoy)/2;
        Hnp(RIndex,RIndex)=Hnp(RIndex,RIndex)-...
              [(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2]'*Px*GaussXW(2,iW)*LE(ElNoy)/2;
        Gnp(RIndex,:)=Gnp(RIndex,:)+[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2]'*Gx*GaussXW(2,iW)*LE(ElNoy)/2;
        VFnp(RIndex)=VFnp(RIndex)+[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2]'*VFx;
    end
 end