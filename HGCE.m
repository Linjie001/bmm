function [Hi,Pi,Gi,VFi]=HGCE(xp,Node,Mat)
% ENCo,xp,GaussXW,GaussPt
%      np=Node(1,1);	En1=Node(2,1);  En2=Node(3,1);    xp=Node(1,2:3)';    
% 	 Den=Mat(1);    
     NU=Mat(3);    SG=Mat(4);
     VFi=zeros(2,1);
     DLE=Node(2,2:3)-Node(1,2:3);
     LEE=sum(DLE.*DLE,2).^0.5;
     l=LEE;
     DLE=DLE./[LEE LEE];
     %NEE=DLE*[0 -1;1 0];
     %H矩阵
     DLR=[xp(1)-Node(1,2) xp(2)-Node(1,3)];
     r1=sum(DLR.*DLR,2).^0.5;
     r2=l-r1;
     m1=r2/l*log(r2/r1)-1;m2=r1/l*log(r2/r1)+1;
     Hi=(1-2*NU)/4/pi/(1-NU)*[0 -m1 0 -m2;m1 0 m2 0];
%      HS2=[r2/l 0 r1/l 0;0 r2/l 0 r1/l]*1/2;
%      H=HS1;%+HS2;
     Pi=(1-2*NU)/4/pi/(1-NU)*[0 -log(r2/r1);log(r2/r1) 0];
     %G矩阵
     n1=r1/2*log(r1)+r2*(r2*log(r2)+r1*log(r1))/l/2-(3*r2+r1)/4;
     n2=r2/2*log(r2)+r1*(r2*log(r2)+r1*log(r1))/l/2-(3*r1+r2)/4;
     GS1=(4*NU-3)*[n1 0 n2 0;0 n1 0 n2];
%      GS2=l/2*[DLE(1)*DLE(1) DLE(1)*DLE(2) DLE(1)*DLE(1) DLE(1)*DLE(2);...
%               DLE(2)*DLE(1) DLE(2)*DLE(2) DLE(2)*DLE(1) DLE(2)*DLE(2)];
     GS2=l/2*[DLE(1)*DLE(1)-(7-8*NU)/2 DLE(1)*DLE(2) DLE(1)*DLE(1)-(7-8*NU)/2 DLE(1)*DLE(2);...
          DLE(2)*DLE(1) DLE(2)*DLE(2)-(7-8*NU)/2 DLE(2)*DLE(1) DLE(2)*DLE(2)-(7-8*NU)/2];
     Gi=1/8/pi/SG/(1-NU)*(GS1+GS2);
	 %%%计算体积力Bi
	%VForce=rl*(2*log(rl)+1)*(Den*AG*([drd1 drd2]*NEE')-...
    %                        1/2/(1-NU)*Den*AG'*[drd1 drd2]'*NEE');
	%VFi=VFi-VForce*GaussXW(2,iG)*LEE/2/8/pi/SG ;                            
    end

