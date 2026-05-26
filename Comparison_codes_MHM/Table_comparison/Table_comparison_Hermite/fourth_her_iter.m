function [t_iter, xc] = fourth_her_iter(n)

epsil = 1e-10;

Nn  = floor(n/2);
ind = n - 2*Nn;

xp = 0;
delta = (1 + ind)/2 * pi / sqrt(2*n + 1);
xn = xp + delta;

xc   = zeros(Nn+ind,1);
deri = zeros(Nn+ind,1);
I    = zeros(Nn+ind,1);

if ind > 0
    u0  = 0.0;
    up1 = 1.0;
    i = 1;
    xc(i) = 0.0;
    deri(i) = 1.0;
    I(i) = 0;
else
    u0  = 1.0;
    up1 = 0.0;
    i = 0;
end

while i < Nn + ind

    i = i + 1;
    err = 1 + epsil;
    I(i) = 0;                      % iteration counter for this root

    while err > epsil

        t = xn - xp;
        delta = t;

        um1 = 0.0;
        um2 = 0.0;

        suma  = u0 + up1*t;
        sumad = up1;

        k = -1;
        errod = 1;

        while ((errod > 1.0e-25 && k < 50) || (k < 5))

            k = k + 1;

            up2 = (xp^2 - 2*n - 1.0)*u0 + 2.0*k*xp*um1 + k*(k-1.0)*um2;

            um2 = um1;
            um1 = u0;
            u0  = up1;
            up1 = up2;

            ntd = up2*t;
            sumad = sumad + ntd;

            if (i > 1+ind)
                errod = abs(ntd/sumad);
            end

            t = t*delta/(k+2.0);

            nt = up2*t;
            suma = suma + nt;

        end

        u0  = suma;
        up1 = sumad;

        hs = u0/up1;
        ws = sqrt(2*n + 1 - xn^2);

        delta = -atan(ws*hs)/ws;

        xp = xn;
        xn = xp + delta;

       I(i) = I(i) + 1;

        err = abs(delta) / max(abs(xn),1);   % safer stopping test

    end

    xc(i) = xn;
   

    % recompute derivative for weights

    t = xn - xp;
    delta = t;

    um1 = 0.0;
    um2 = 0.0;

    suma  = u0 + up1*t;
    sumad = up1;

    k = -1;
    errod = 1.0;

    while ((errod > 1.0e-25 && k < 50) || (k < 5))

        k = k + 1;

        up2 = (xp^2 - 2*n - 1.0)*u0 + 2*k*xp*um1 + k*(k-1)*um2;

        um2 = um1;
        um1 = u0;
        u0  = up1;
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

    u0  = suma;
    up1 = sumad;

    deri(i) = sumad;

    delta = pi / sqrt(2*n + 1 - xn^2);

    xp = xn;
    xn = xp + delta;

end


% -------- Compute weights --------

finx = 20.0;
sumap = 0.0;

k = 1;

while (k <= Nn+ind && xc(k) < finx)

    xp2 = xc(k)^2;

    sumap = sumap + 2*xp2*exp(-xp2)/deri(k)^2;

    k = k + 1;

end


w   = zeros(Nn+ind,1);
wns = zeros(Nn+ind,1);

for k = 1:Nn+ind

    w(k)   = 0.5*sqrt(pi)/sumap/deri(k)^2;

    wns(k) = w(k)*exp(-xc(k)^2);

end


% total iteration count
t_iter = sum(I);

end