function [x, w] = legendre_fourth_test(n)
  [x1,w1]=fourth_legendre(n);
  x2 = x1;
  w2 = w1;
  s = 0;

  if (~isempty(x1))&&(~isempty(x2))
     if abs(x2(1)+x1(1))<3.0e-16
        s = 1; 
     end
  end   
 x = [-flipud(x2(1+s:end)); x1];
 w = [flipud(w2(1+s:end)); w1];
 le = length(x);
 w1 = w(1:le);
 w1 = w1/w1(1);
 la = 2/sum(w1);
 w = la*w1;

function [x,w] = fourth_legendre(n)
x = [];
w = [];
Nn = floor(n/2);
ind = n-2*Nn;
x = zeros(Nn+ind, 1);
w = zeros(Nn+ind, 1);
L = 2*n+1;
nor = 1.e-130;
ind = n-2*floor(n/2);
   if ind>0
       y0=0;y1=nor; i=1; x(1)=0; w(1)=1;
   else
       y0=nor;y1=0; i=0;
   end
   xi= tanh(((ind+1)*pi)/sqrt(L^2-1));
   xt=0;
   t=xi-xt;
   [y0, y1] = seriestay(n, xt, t, y0, y1);
   [xt, wt, y0, y1] = fourth_ser(n, t, y0, y1);
   i = i + 1;
   x(i) = xt;
   w(i) = wt;
   sqde=sqrt(0.25*((L*L-1)*(1-xt*xt)));
   de=-tanh(pi/sqde);
   xi=(xt-de)/(1-xt*de);
   t=xi-xt;
   [y0,y1]=seriestay(n,xt,t,y0,y1);
   while i < Nn+ind
       [xt,wt,y0,y1] = fourth_ser(n,xi,y0,y1);
       i = i+1;
       x(i) = xt;
       w(i) = wt;
       sqde = sqrt(0.25*(L*L-1)*(1-xt*xt));
       de=-tanh(pi/sqde);
       xi=(xt-de)/(1-xt*de);
       t=xi-xt;
       [y0,y1]=seriestay(n,xt,t,y0,y1);
   end

 
function [x,w,y0,y1]=fourth_ser(n,x0,y0,y1)
% fixed point on z and using Taylor for (1-x)^(a+1)/2*(1+x)^(b+1)/2*P
epsil=10^(-15);
erro=1;
L=2*n+1;
% h = y0/(y1*(1-x0*x0)+x0*y0);
% while ~h>0
%     sq = 0.25*(L*L-1)*(1-x0*x0);
%     sqde = sqrt(sq);
%     xp = x0;
%     argu=atan(sqde*h)/sqde-pi/sqde;
%     de = tanh(argu);
%     x0 = (x0-de)/(1-x0*de);
%     t = x0 - xp;
%     [y0,y1]=seriestay(n,x0,t,y0,y1);
%     h = y0/(y1*(1-x0*x0)+x0*y0);
% end

while erro>epsil
   sqde=sqrt(0.25*(L*L-1)*(1-x0*x0));
   h=y0/(y1*(1-x0*x0)+x0*y0);
   argu=atan(sqde*h)/sqde;    
   if ~h>0
      argu=argu-pi/sqde; 
   end    
      de=tanh(argu);
   x=(x0-de)/(1-x0*de);   
   erro=abs(1-x/x0);
   t=x-x0;
   [y0,y1]=seriestay(n,x0,t,y0,y1);
   x0=x;
end
w=1/y1/y1;




function [y0,yp1]=seriestay(n,xp,t,y0,yp1)
% TAYLOR EN X para (1-x)^((a+1)/2)*(1+x)*((b+1)/2)P(n,a,b,x)
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
     ym2=ym1;ym1=y0;y0=yp1;yp1=yp2;  
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

