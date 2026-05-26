
function [x,x0,y0,y1, i]=fixedzser_halley_1(n,x0,y0,y1)
% FOM on z and using Taylor for (1-x)^(1)/2*(1+x)^(1)/2*P 

    L=2*n+1;
    epsil=10^(-15);
    sq=0.25*((L*L-1)*(1-x0*x0));
   
    h=y0/(y1*(1-x0*x0)+x0*y0);
   
    argu=(2*h)/(2+sq*h^2);       
   de=tanh(argu);
   x=(x0-de)/(1-x0*de);   
   erro=abs(1-x/x0);
   i = 1;
    while erro>epsil
        t=x-x0;
        [y0,y1]=seriestay_1(n,x0,t,y0,y1);
        x0=x;
        h=y0/(y1*(1-x0*x0)+x0*y0);   
        argu=(2*h)/(2+sq*h^2);       
        de=tanh(argu);
        x=(x0-de)/(1-x0*de);   
        erro=abs(1-x/x0);
        i = i+1;
    end
end

