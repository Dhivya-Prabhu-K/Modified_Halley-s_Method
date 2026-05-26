
function [y0,yp1]=seriestay_1(n,xp,t,y0,yp1)
% TAYLOR EN X para (1-x)^((a+1)/2)*(1+x)^((b+1)/2)*P(n,a,b,x)
% for the case of legendre a=b=0
 epso=1.0e-19;
 ym1=0;ym2=0; 
 delta=t;
 y0=y0*1.e-30;
 yp1=yp1*1.e-30;
 suma=y0+yp1*t;
 sumad=yp1;
 j=-1;
 errod=1; 
 xp2=xp*xp;
 
 L=2*n+1;L2=L*L;
 while (errod>epso)&&(j<100)
     j=j+1;
     c1=4*(j+2)*(j+1)*(1-xp2)^2;
     c2=-16*(j+1)*j*(1-xp2)*xp;
     c3=0.5*j*(j-1)*(48*xp2-16)+(L2-1)*(1-xp2)+(2)*(1-xp)+(2)*(1+xp);
     c4=16*(j-1)*(j-2)*xp-2*(L2-1)*xp;
     c5=4*(j-2)*(j-3)+(1-L2);
     yp2= -1/c1*(c2*yp1+c3*y0+c4*ym1+c5*ym2);
     ym2=ym1;   ym1=y0;     y0=yp1;     yp1=yp2;  
     ntd=(j+2)*yp2*t;
     sumad=sumad+ntd;
     t=t*delta;
     nt=yp2*t;
     suma=suma+nt;
     if (j>20) && (y0*yp1~=0)
        errod=max(abs((ntd/sumad)),abs((nt/suma)));
     end
 end
 y0=suma;yp1=sumad;
 y0=y0*1.e30;
 yp1=yp1*1.e30;
end
