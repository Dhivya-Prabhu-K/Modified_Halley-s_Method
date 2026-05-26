function [x, q, I]=fsweep_fourth_1(n)

ind=n-2*floor(n/2);
x = zeros(floor(n/2)+ind, 1);
I = zeros(floor(n/2)+ind, 1);
L=2*n+1;
nor=1.e-130;
ind=n-2*floor(n/2);
   if ind>0
       y0=0;y1=nor; i=1; x(1)=0;
   else
       y0=nor;y1=0; i=0;
   end

      xi= tanh(((ind+1)*pi)/sqrt(L^2-1));
      xt=0;
      t=xi-xt;
      [y0, y1] = seriestay_1(n, xt, t, y0, y1);
      [xt, x0, y0, y1, k] = fixedzser_fourth_1(n, t, y0, y1);
        i = i + 1;
        x(i) = xt;
        I(i) = k;
   sqde=sqrt(0.25*((L*L-1)*(1-xt*xt)));
   de=-tanh(pi/sqde);
   xi=(xt-de)/(1-xt*de);
   t=xi-x0;
   [y0,y1]=seriestay_1(n,x0,t,y0,y1);
while i<floor(n/2)+ind
   [xt,x0,y0,y1, k]=fixedzser_fourth_1(n,xi,y0,y1);
   i=i+1;
   x(i)=xt;   
   I(i)=k;
   sqde=sqrt(0.25*((L*L-1)*(1-xt*xt)));
   de=-tanh(pi/sqde);
   xi=(xt-de)/(1-xt*de);
   t=xi-x0;
   [y0,y1]=seriestay_1(n,x0,t,y0,y1);
end
q = sum(I);
end