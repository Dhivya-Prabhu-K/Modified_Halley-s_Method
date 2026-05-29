function [x, w] = hermite_asy(n)
    % double precision
    x0 = hermite_initialguess(n);     % initial guesses
    t0 = x0 ./ sqrt(2*n+1);
    theta = acos(t0);                 % θ-variable
    val = x0;

    % Newton iterations in θ-space
    for k = 1:20
        [val0, val1] = hermpoly_asy_airy(n, theta);
        dtheta = -val0 ./ (sqrt(2) * sqrt(2*n+1) .* val1 .* sin(theta));
        theta = theta - dtheta;
        if norm(dtheta, inf) < sqrt(eps) / 10
            break;
        end
    end

    t0 = cos(theta);
    x = sqrt(2*n+1) * t0;             % back to x-variable
    w = x .* val0 + sqrt(2) * val1;
    w = 1 ./ (w.^2);                  % quadrature weights

    %Symmetrize nodes and weights
    if mod(n,2) == 1
        w_full = [flip(w(:)); w(2:end)];
        x_full = [-flip(x(:)); x(2:end)];
    else
        w_full = [flip(w(:)); w(:)];
        x_full = [-flip(x(:)); x(:)];
    end

    w_full = w_full * sqrt(pi) / sum(exp(-x_full.^2).*w_full);
    x = x_full;
    w = w_full;
    % if mod(n,2) == 1
    %     s=2*sum(exp(-x.^2).*w)-w(1);
    % else
    %     s=2*sum(exp(-x.^2).*w);
    % end
    % w=w*sqrt(pi)/s;
end


function [val, dval] = hermpoly_asy_airy(n, theta)
    % Evaluate Hermite polynomial asymptotics in theta-space (Airy-based)
    
    musq = 2*n+1;
    cosT = cos(theta);
    sinT = sin(theta);
    sin2T = 2*cosT.*sinT;
    
    eta = 0.5*theta - 0.25*sin2T;
    chi = -(3*eta/2).^(2/3);
    phi = (-chi ./ sinT.^2).^(1/4);
    C = 2*sqrt(pi)*musq^(1/6).*phi;
    z = musq^(2/3).*chi;
    % Airy functions (MATLAB built-ins)

    Ai = real(airy(0, z));
    Aip = real(airy(1, z));
    % u-polynomials
    u0 = 1;
    u1 = (cosT.^3 - 6*cosT)/24;
    u2 = (-9*cosT.^4 + 249*cosT.^2 + 145)/1152;
    u3 = (-4042*cosT.^9 + 18189*cosT.^7 - 28287*cosT.^5 ...
        - 151995*cosT.^3 - 259290*cosT)/414720;

    % coefficients (a,b)
    a0=1; b0=1;
    a1=15/144; b1=-7/5*a1;
    a2=5*7*9*11/2/144^2; b2=-13/11*a2;
    a3=7*9*11*13*15*17/6/144^3; b3=-19/17*a3;

    % Hermite polynomial value
    val = Ai;
    B0 = -(a0*u1.*phi.^6 + a1*u0) ./ chi.^2;
    val = val + B0.*Aip/musq^(4/3);

    A1 = (b0*u2.*phi.^12 + b1*u1.*phi.^6 + b2*u0) ./ chi.^3;
    val = val + A1.*Ai/musq^2;

    B1 = -(u3.*phi.^18 + a1*u2.*phi.^12 + a2*u1.*phi.^6 + a3*u0) ./ chi.^5;
    val = val + B1.*Aip/musq^(4/3+2);

    val = C.*val;

    % Derivative part

    C = sqrt(2*pi)*musq^(1/3)./phi;

    % v-polynomials
    v0 = 1;
    v1 = (cosT.^3 + 6*cosT)/24;
    v2 = (15*cosT.^4 - 327*cosT.^2 - 143)/1152;
    v3 = (259290*cosT + 238425*cosT.^3 - 36387*cosT.^5 ...
        + 18189*cosT.^7 - 4042*cosT.^9)/414720;

    % derivative
    C0 = -(b0*phi.^6.*v1 + b1*v0)./chi;
    dval = C0.*Ai/musq^(2/3);

    D0 = a0*v0;
    dval = dval + D0*Aip;

    C1 = -(v3.*phi.^18 + b1*v2.*phi.^12 + b2*v1.*phi.^6 + b3*v0)./chi.^4;
    dval = dval + C1.*Ai/musq^(2/3+2);

    D1 = (a0*v2.*phi.^12 + a1*v1.*phi.^6 + a2*v0)./chi.^3;
    dval = dval + D1.*Aip/musq^2;

    dval = C.*dval;
end


function x_init = hermite_initialguess(n)
    % Initial guesses for Hermite polynomial zeros (Gatteschi + Tricomi)
    if n < 20
        error('n must be >= 20 for asymptotic initial guesses');
    end

    if mod(n,2)==1
        m = (n-1)/2;
        a = 0.5;
    else
        m = n/2;
        a = -0.5;
    end
    nu = 4*m + 2*a + 2;

    T = @(t) t.^(2/3).*(1 + 5/48*t.^(-2) - 5/36*t.^(-4) + ...
        (77125/82944)*t.^(-6) - 108056875/6967296*t.^(-8) + ...
        162375596875/334430208*t.^(-10));

    airyrts = -T(3*pi/8*(4*(1:m)-1));


    % TODO: replace first 10 roots with precomputed AIRY_ROOTS if available
    airy_roots = [-2.338107410459767, -4.087949444130970, -5.520559828095551, ...
                       -6.786708090071759, -7.944133587120853, -9.022650853340980, ...
                       -10.04017434155809, -11.00852430373326, -11.93601556323626, ...
                       -12.82877675286576];
    airyrts(1:min(10, m)) = airy_roots(1:min(10, m));

    x_init_airy = sqrt(abs(nu + (2^(2/3))*airyrts*nu^(1/3) + ...
        (1/5*2^(4/3))*airyrts.^2*nu^(-1/3) + ...
        (11/35-a^2-12/175)*airyrts.^3/nu));

    x_init_airy = flip(x_init_airy(:));

    % Tricomi-type initial guesses near 0
    Tnk0 = pi/2*ones(m,1);
    rhs = ((4*m+3) - 4*(1:m))/nu * pi;
    for k = 1:7
        val = Tnk0 - sin(Tnk0) - rhs(:);
        dval = 1 - cos(Tnk0);
        Tnk0 = Tnk0 - val./dval;
    end
    tnk = cos(Tnk0/2).^2;
    x_init_sin = sqrt(nu*tnk - ((tnk+1/4)./(tnk-1).^2 + (3*a^2-1))/(3*nu));

    % Patch together
    p = 0.4985+eps;
    idx1 = floor(p*n);
    idx2 = ceil(p*n);
    x_init = [x_init_sin(1:idx1); x_init_airy(idx2:end)];

    if mod(n,2)==1
        x_init = [0; x_init];
        x_init = x_init(1:m+1);
    else
        x_init = x_init(1:m);
    end
end
