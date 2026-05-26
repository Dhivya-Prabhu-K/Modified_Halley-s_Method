function [t_iter, xc] = third_her_iter(n)
epsil = 1e-10;
Nn = floor(n/2);
ind = n - 2*Nn;

     xp = 0;
    delta = pi / (2*sqrt(2*(n+1)));
    xn = xp + delta;
    xc   = zeros(Nn+ind,1);
    deri = zeros(Nn+ind,1);
    I = zeros(Nn+ind,1);
    if ind > 0
        u0 = 0;
        up1 = 1;
        i = 1;
        xc(i) = 0;
        deri(i) = 1;
        I(i) = 0;
    else
        u0 = 1; up1 = 0; i = 0;
        
    end
    while i < 3
        i = i + 1;
        err = 1 + epsil;
        I(i) = 0;
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
                up2 = (xp^2 - 2*n - 1.0)*u0 + 2.0*k*xp*um1 + k*(k-1.0)*um2;
                um2 = um1;
                um1 = u0;
                u0 = up1;
                up1 = up2;

                ntd = up2*t;
                sumad = sumad + ntd;
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
            I(i) = I(i) + 1;
            err = abs(delta / xn);
        end
        xc(i) = xn;
        I(i) = j;
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
        I(i) = 0;
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
            I(i) = I(i) + 1;
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
t_iter = sum(I);
end
