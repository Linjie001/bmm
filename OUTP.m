%function [Gi Hi AD AS FV SV]=OUTP(GaussXW,Node([np En1 En2]),Mat(MatiA,:),xp,AG)
function [Gi,Hi,AD,AS,FV,SV]=OUTP(Node,Mat,xp,AG,GaussPt)   
   
[x, w] = GaussLegendre(GaussPt);
GaussXW = [x w]';
	Den=Mat(1);    NU=Mat(3);    SG=Mat(4);
    
    En1=Node(1,1);En2=Node(2,1);
                AD=zeros(3,4);
                AS=zeros(3,4);
                Gi=zeros(2,4);
                Hi=zeros(2,4);
                SV=zeros(3,1);
                SV11=0;
                SV12=0;
                SV22=0;
                FV=0;
                 DL=Node(2,2:3)-Node(1,2:3);
                 LE=sum(DL.*DL,2).^0.5;
                 DL=DL./[LE LE];
                 NE=DL*[0 -1;1 0];
               for iG=1:size(GaussXW,2)
                    x=GaussXW(1,iG);
                    xl=([1-x x+1]/2*Node([1 2],2:3)).';
                    r=xl-xp;
                    rl=sqrt(r.'*r)  ;        
                    drd1=r(1)/rl;
                    drd2=r(2)/rl;
                    drdn=[drd1 drd2]*NE';
                    NShape=[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2];
        %             计算AU
                    u11=(4*NU-3)*log(rl)+drd1^2-(7-8*NU)/2;
                    u12=drd1*drd2;
                    u21=u12;
                    u22=(4*NU-3)*log(rl)+drd2^2-(7-8*NU)/2;
                    U=[u11 u12;u21 u22]*1/8/pi/SG/(1-NU)*LE/2;  
                    Gi=Gi+U*NShape*GaussXW(2,iG); 
        %              计算AG
                    p11=drdn*((1-2*NU)+2*drd1^2);
                    p12=drdn*(2*drd1*drd2)-(1-2*NU)*(drd1*NE(2)-drd2*NE(1));
                    p21=drdn*(2*drd1*drd2)+(1-2*NU)*(drd1*NE(2)-drd2*NE(1));
                    p22=drdn*((1-2*NU)+2*drd2^2); 
                    P=-1/4/pi/(1-NU)/rl*LE/2*[p11 p12;p21 p22];
                    Hi=Hi+P*NShape*GaussXW(2,iG);
        %            计算U（位移）中的FVjj
                    VForce=rl*(2*log(rl)+1)*(Den*AG*([drd1 drd2]*NE')-...
                           1/2/(1-NU)*Den*AG'*[drd1 drd2]'*NE');
                    FV=FV-VForce*GaussXW(2,iG)*1/8/pi/SG*LE/2;    
         %             计算AD
                    D111=(1-2*NU)*drd1+2*drd1^3;
                    D112=(1-2*NU)*drd2+2*drd1^2*drd2;
                    D122=(1-2*NU)*(-drd1)+2*drd1*drd2^2;
                    D211=(1-2*NU)*(-drd2)+2*drd2*drd1^2;
                    D212=(1-2*NU)*drd1+2*drd2^2*drd1;
                    D222=(1-2*NU)*drd2+2*drd2^3;        
                    D=1/4/pi/(1-NU)/rl*LE/2*[D111 D112 D122;D211 D212 D222]';  
                    AD=AD+D*NShape*GaussXW(2,iG);
        %              计算AS
%                     a1=2*drdn*((1-2*NU)*drd1+NU*(2*drd1));
%                     a2=2*NU*(2*drd1^2*NE(ElNo,1))+(1-2*NU)*(2*drd1^2*NE(ElNo,1)+2*NE(ElNo,1))-(1-4*NU)*NE(ElNo,1);
                    S111=2*drdn*((1-2*NU)*drd1+NU*(2*drd1)-4*drd1^3)+2*NU*(2*drd1^2*NE(1))+...
                        (1-2*NU)*(2*drd1^2*NE(1)+2*NE(1))-(1-4*NU)*NE(1);
                    S112=2*drdn*(NU*drd2-4*drd1^2*drd2)+2*NU*(drd1*drd2*NE(1)+drd1^2*NE(2))+...
                        (1-2*NU)*(2*drd1*drd2*NE(1)+NE(2));
                    S122=2*drdn*((1-2*NU)*drd1-4*drd1*drd2^2)+2*NU*(2*drd1*drd2*NE(2))+...
                        (1-2*NU)*(2*drd2^2*NE(1))-(1-4*NU)*NE(1);

                    S211=2*drdn*((1-2*NU)*drd2-4*drd1^2*drd2)+2*NU*(2*drd1*drd2*NE(1))+...
                        (1-2*NU)*(2*drd1^2*NE(2))-(1-4*NU)*NE(2);
                    S212=2*drdn*(NU*drd1-4*drd2^2*drd1)+2*NU*(drd2^2*NE(1)+drd2*drd1*NE(2))+...
                        (1-2*NU)*(2*drd1*drd2*NE(2)+NE(1));
                    S222=2*drdn*((1-2*NU)*drd2+NU*(2*drd2)-4*drd2^3)+2*NU*(2*drd2^2*NE(2))+...
                        (1-2*NU)*(2*drd2^2*NE(2)+2*NE(2))-(1-4*NU)*NE(2);      
                    S=[S111 S112 S122;S211 S212 S222]'*SG/2/pi/(1-NU)/rl^2*LE/2;
                    AS=AS+S*NShape*GaussXW(2,iG);
                    
        %              计算P（应力）中的SVjj
                    SV11=SV11+1/8/pi*(LE/2)*GaussXW(2,iG)*(1/(1-NU)*(NU*(2*NE*[drd1 drd2]'*Den*AG'*[drd1 drd2]'+...
                        (1-2*log(1/rl))*Den*AG'*NE')-Den*AG'*[drd1 drd2]'*(2*NE(1)*drd1)));
                    SV12=SV12+1/8/pi*(LE/2)*GaussXW(2,iG)*(2*NE*[drd1 drd2]'*(Den*AG(2)*drd1)+...
                        1/(1-NU)*(-Den*AG'*[drd1 drd2]'*(NE*[drd2 drd1]')+...
                        (1-2*NU)/2*(1-2*log(1/rl))*(Den*AG(2)*NE(1))));
                    SV22=SV22+1/8/pi*(LE/2)*GaussXW(2,iG)*(2*NE*[drd1 drd2]'*(2*Den*AG(2)*drd2)+1/(1-NU)*(NU*((2*NE*[drd1 drd2]')*Den*AG'*[drd1 drd2]'+...
                        (1-2*log(1/rl))*(Den*AG'*NE'))-Den*AG'*[drd1 drd2]'*(2*NE(2)*drd2)+(1-2*NU)/2*(1-2*log(1/rl))*(2*Den*AG(2)*NE(2))));
                    SV=[SV11 SV12 SV22]';
               end