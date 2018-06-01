function [Hi,Pi,Gi,VFi]=HGFE(xp,GaussXW,Node,Mat,AG)
% ENCo,xp,GaussXW,GaussPt
%     En1=Node(1,1);  En2=Node(2,1);    
	Den=Mat(1);    NU=Mat(3);    SG=Mat(4);
    Hi=zeros(2,4);Pi=zeros(2);Gi=zeros(2,4);VFi=zeros(2,1);
    DLE=Node(2,2:3)-Node(1,2:3);
    LEE=sum(DLE.*DLE,2).^0.5;
    DLE=DLE./[LEE LEE];
    NEE=DLE*[0 -1;1 0];
    
	for iG=1:size(GaussXW,2)
        x=GaussXW(1,iG);
        xl=([1-x x+1]/2*Node(1:2,2:3)).'; 
        r=xl-xp;
        rl=sqrt(r.'*r);            
        drd1=r(1)/rl;
        drd2=r(2)/rl;
        drdn=[drd1 drd2]*NEE';
        NShape=[(1-x)/2 0 (1+x)/2 0;0 (1-x)/2 0 (1+x)/2];
        %求H矩阵
        p11=drdn*((1-2*NU)+2*drd1^2);
        p12=drdn*(2*drd1*drd2)-(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
        p21=drdn*(2*drd1*drd2)+(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
        p22=drdn*((1-2*NU)+2*drd2^2); 
        P=-1/4/pi/(1-NU)/rl*[p11 p12;p21 p22]*LEE/2;
        FPN=P*NShape;
        Hi=Hi+FPN*GaussXW(2,iG);
        Pi=Pi+P*GaussXW(2,iG);
        %求G矩阵
        u11=(4*NU-3)*log(rl)+drd1^2-(7-8*NU)/2 ;     %;               
        u12=drd1*drd2;
        u21=u12;
        u22=(4*NU-3)*log(rl)+drd2^2-(7-8*NU)/2 ;     %;  
        U=1/8/pi/SG/(1-NU)*LEE/2*[u11 u12;u21 u22];
        Gi=Gi+U*NShape*GaussXW(2,iG); 
	    %%%计算体积力Bi
	    VForce=rl*(2*log(rl)+1)*(Den*AG*([drd1 drd2]*NEE')-...
                            1/2/(1-NU)*Den*AG'*[drd1 drd2]'*NEE');
	    VFi=VFi-VForce*GaussXW(2,iG)*LEE/2/8/pi/SG ;                            
    end