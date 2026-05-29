function [x_value, y_value, x1_value, y1_value] = error_third_her(n)
% double precision
epsil = 1e-10;
Nn = floor(n/2);
ind = n - 2*Nn;

     xp = 0;
    delta = pi / (2*sqrt(2*(n+1)));
    xn = xp + delta;
     xc   = zeros(Nn+ind,1);
    deri = zeros(Nn+ind,1);
    if ind > 0
        u0 = 0; up1 = 1; i = 1;
        xc(i) = 0.0;
        deri(i) = 1.0;
    else
        u0 = 1; up1 = 0; i = 0;
        
    end
    while i < 3
        i = i + 1;
        err = 1 + epsil;
        while err > epsil
            t = xn - xp;
            delta = t;
            um1 = 0; um2 = 0;
            suma = u0 + up1 * t;
            sumad = up1;
            k = -1;
            errod = 1;
            
           
            while (errod >1e-25 && k < 50) || k < 5
                k = k + 1;
                up2 = (xp^2 - 2*n - 1)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
                um2 = um1; um1 = u0; u0 = up1; up1 = up2;
                ntd = up2 * t; sumad = sumad + ntd;
                if i > 1 + ind
                    errod = abs(ntd / sumad);
                end
                t = t * delta / (k + 2);
                nt = up2 * t;
                suma = suma + nt;
            end
            u0 = suma; up1 = sumad;
            hs = -sqrt(2*(n + 1))*(u0 / (xn*u0 - up1));
            ws = -xn/sqrt(2*(n+1)) ;
            delta = -(2*hs) / (sqrt(2*(n + 1))*(2 + hs^2 - 2*ws*hs));
            xp = xn;
            xn = xp + delta;
            err = abs(delta / xn);
        end
        xc(i) = xn;

          % recompute derivative for weights
        t = xn - xp;
        delta = t;
        um1 = 0.0;
        um2 = 0.0;
        suma = u0 + up1*t;
        sumad = up1;
        k = -1;
        errod = 1.0;
        while ((errod > 1.0e-25 && k < 50) || (k < 5))
            k = k + 1;
            up2 = (xp^2 - 2*n - 1.0)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
            um2 = um1;
            um1 = u0;
            u0 = up1;
            up1 = up2;
            ntd = up2*t;
            sumad = sumad + ntd;
            if (i > 1+ind)
                errod = abs(ntd/sumad);
            end
            t = t*delta/(k+2);
            nt = up2*t;
            suma = suma + nt;
        end
        u0 = suma;
        up1 = sumad;
        deri(i) = sumad;
        xp = xn;
        xn = xp + pi/(2*sqrt(2*(n+1)));
    end
    while i >= 3 && i < Nn + ind
        i = i + 1;
        err = 1 + epsil;
        while err > epsil
            t = xn - xp;
            delta = t;
            um1 = 0; um2 = 0;
            suma = u0 + up1 * t;
            sumad = up1;
            k = -1;
            errod = 1;
      
            while (errod > 1e-25 && k < 50) || k < 5
                k = k + 1;
                up2 = (xp^2 - 2*n - 1)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
                um2 = um1; um1 = u0; u0 = up1; up1 = up2;
                ntd = up2 * t; sumad = sumad + ntd;
                if i > 1 + ind
                    errod = abs(ntd / sumad);
                end
                t = t * delta / (k + 2);
                nt = up2 * t;
                suma = suma + nt;
            end
            u0 = suma; up1 = sumad;
            hs = -sqrt(2*(n + 1))*(u0 / (xn*u0 - up1));
            ws = -xn/sqrt(2*(n+1)) ;
            delta = -(2*hs) / (sqrt(2*(n + 1))*(2 + hs^2 - 2*ws*hs));
            xp = xn;
            xn = xp + delta;
            err = abs(delta / xn);
        end
        xc(i) = xn;
          % recompute derivative for weights
        t = xn - xp;
        delta = t;
        um1 = 0.0;
        um2 = 0.0;
        suma = u0 + up1*t;
        sumad = up1;
        k = -1;
        errod = 1.0;
        while ((errod > 1.0e-25 && k < 50) || (k < 5))
            k = k + 1;
            up2 = (xp^2 - 2*n - 1.0)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
            um2 = um1;
            um1 = u0;
            u0 = up1;
            up1 = up2;
            ntd = up2*t;
            sumad = sumad + ntd;
            if (i > 1+ind)
                errod = abs(ntd/sumad);
            end
            t = t*delta/(k+2);
            nt = up2*t;
            suma = suma + nt;
        end
        u0 = suma;
        up1 = sumad;
        deri(i) = sumad;
        xp = xn;
        xn = xp + (xc(i-1) - xc(i-2));
    end


    % Compute weights
    finx = 20.0;
    sumap = 0.0;
    i = 0;
    while (i < Nn+ind && xc(max(i,1)) < finx)
        i = i + 1;
        xp2 = xc(i)^2;
        sumap = sumap + 2*xp2*exp(-xp2)/deri(i)^2;
    end

    w   = zeros(Nn+ind,1);
    wns = zeros(Nn+ind,1);
    for i = 1:Nn+ind
        w(i)   = 0.5*sqrt(pi)/sumap/deri(i)^2;
        wns(i) = w(i)*exp(-xc(i)^2);
    end
    x16=xc;
    w16=w;

    
    % extended precision third_her method

digits(32)
n = vpa(n);
epsil = vpa((10)^(-11));
Nn = floor(n/2);
ind = n - 2*Nn;

xp = vpa(0);
    delta = vpa(pi) / vpa((2*sqrt(2*(n+1))));
    xn = xp + delta;
    % Preallocation
    xc   = sym(zeros(double(Nn+ind),1));
    deri = sym(zeros(double(Nn+ind),1));
    if ind > 0
        u0 = vpa(0);
        up1 = vpa(1);
        i = 1;
        xc(i) = vpa(0);
        deri(i) = vpa(1);
    else
        u0 = vpa(1);
        up1 = vpa(0);
        i = 0;
        
    end
    while i < 3
        i = i + 1;
        err = 1 + epsil;
        while err > epsil
            t = xn - xp;
            delta = t;
        um1 = vpa(0); um2 = vpa(0);
        suma = u0 + up1 * t;
        sumad = up1;
        k = -1;
        errod = vpa(1);
        tol = vpa(1e-80);
        max_k = 100;
           
            while (errod >tol && k < max_k) || k < 5
                k = k + 1;
                up2 = (xp^2 - 2*n - 1)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
                um2 = um1; um1 = u0; u0 = up1; up1 = up2;
                ntd = up2 * t; sumad = sumad + ntd;
                if i > 1 + ind
                    errod = abs(ntd / sumad);
                end
                t = t * delta / (k + 2);
                nt = up2 * t;
                suma = suma + nt;
            end
            u0 = suma; up1 = sumad;
            hs = -vpa(sqrt(2*(n + 1)))*(u0 / (xn*u0 - up1));
            ws = -xn/vpa(sqrt(2*(n+1))) ;
            delta = -(2*hs) / vpa(sqrt(2*(n + 1))*(2 + hs^2 - 2*ws*hs));
            xp = xn;
            xn = xp + delta;
            err = abs(delta / xn);
        end
        xc(i) = xn;

          % recompute derivative for weights
        t = xn - xp;
        delta = t;
        um1 = vpa(0); um2 = vpa(0);
        suma = u0 + up1 * t;
        sumad = up1;
        k = -1;
        errod = vpa(1);
        tol = vpa(1e-80);
        max_k = 100;
        while ((errod > tol && k < max_k) || (k < 5))
            k = k + 1;
            up2 = (xp^2 - 2*n - 1.0)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
            um2 = um1;
            um1 = u0;
            u0 = up1;
            up1 = up2;
            ntd = up2*t;
            sumad = sumad + ntd;
            if (i > 1+ind)
                errod = abs(ntd/sumad);
            end
            t = t*delta/(k+2);
            nt = up2*t;
            suma = suma + nt;
        end
        u0 = suma;
        up1 = sumad;
        deri(i) = sumad;
        xp = xn;
        xn = xp + vpa(pi)/vpa(2*sqrt(2*(n+1)));
    end
    while i >= 3 && i < Nn + ind
        i = i + 1;
        err = 1 + epsil;
        while err > epsil
            t = xn - xp;
            delta = t;
                  um1 = vpa(0); um2 = vpa(0);
        suma = u0 + up1 * t;
        sumad = up1;
        k = -1;
        errod = vpa(1);
        tol = vpa(1e-80);
        max_k = 100;
            while ((errod > tol && k < max_k) || (k < 5))
                k = k + 1;
                up2 = (xp^2 - 2*n - 1)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
                um2 = um1; um1 = u0; u0 = up1; up1 = up2;
                ntd = up2 * t; sumad = sumad + ntd;
                if i > 1 + ind
                    errod = abs(ntd / sumad);
                end
                t = t * delta / (k + 2);
                nt = up2 * t;
                suma = suma + nt;
            end
            u0 = suma; up1 = sumad;
            hs = -vpa(sqrt(2*(n + 1)))*(u0 / (xn*u0 - up1));
            ws = -xn/vpa(sqrt(2*(n+1))) ;
            delta = -(2*hs) / vpa(sqrt(2*(n + 1))*(2 + hs^2 - 2*ws*hs));
            xp = xn;
            xn = xp + delta;
            err = abs(delta / xn);
        end
        xc(i) = xn;
          % recompute derivative for weights
        t = xn - xp;
        delta = t;
        um1 = vpa(0); um2 = vpa(0);
        suma = u0 + up1 * t;
        sumad = up1;
        k = -1;
        errod = vpa(1);
        tol = vpa(1e-80);
        max_k = 100;
        while ((errod > tol && k < max_k) || (k < 5))
            k = k + 1;
            up2 = (xp^2 - 2*n - 1.0)*u0 + 2*k*xp*um1 + k*(k-1)*um2;
            um2 = um1;
            um1 = u0;
            u0 = up1;
            up1 = up2;
            ntd = up2*t;
            sumad = sumad + ntd;
            if (i > 1+ind)
                errod = abs(ntd/sumad);
            end
            t = t*delta/(k+2);
            nt = up2*t;
            suma = suma + nt;
        end
        u0 = suma;
        up1 = sumad;
        deri(i) = sumad;
        xp = xn;
        xn = xp + (xc(i-1) - xc(i-2));
    end


    % Compute weights
    finx = 20.0;
    sumap = 0.0;
    i = 0;
    while (i < Nn+ind && xc(max(i,1)) < finx)
        i = i + 1;
        xp2 = xc(i)^2;
        sumap = sumap + 2*xp2*exp(-xp2)/deri(i)^2;
    end

    w   = sym(zeros(double(Nn+ind),1));
    wns = sym(zeros(double(Nn+ind),1));
    for i = 1:Nn+ind
        w(i)   = 0.5*vpa(sqrt(pi))/sumap/deri(i)^2;
        wns(i) = w(i)*vpa(exp(-xc(i)^2));
    end

    x32=xc;
    w32=w;
    % computing the error values
    error = zeros(1, double(Nn + ind));
for j = 1:double(Nn + ind)
    error(j) = abs(1 - x16(j) / (x32(j)));
end
x_value = x16;  % your x values
for i=1:length(xc)
    if error(i)>0
     y_value(i)=log(error(i));
    else
      y_value(i)=-50;
    end
end 
    error1 = zeros(1, double(Nn + ind));
for j = 1:double(Nn + ind)
    error1(j) = abs(1 - w16(j) / (w32(j)));
end
x1_value = w16;  % your x values
for i=1:length(xc)
    if error1(i)>0
     y1_value(i)=log(error1(i));
    else
      y1_value(i)=-50;
    end
end
end
