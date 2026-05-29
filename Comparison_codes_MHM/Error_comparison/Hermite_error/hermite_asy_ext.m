function [x, w] = hermite_asy_ext(n)
    % Keep digits for vpa calculations
    digits(32);

    % Keep a numeric copy for sizes and indexing
    n_num = double(n);

    % Use initial guesses computed in numeric (fast) then convert once
    x0_num = hermite_initialguess_num(n_num);     % numeric initial guesses (double)
    x0 = vpa(x0_num);                              % convert once to vpa

    % Prepare theta in vpa
    t0 = x0 ./ vpa(sqrt(2*n_num+1));
    theta = acos(t0);                 % theta in vpa
    val = x0;

    % Newton iterations in theta-space (vpa arithmetic)
    maxit = 15;
    for k = 1:maxit
        [val0, val1] = hermpoly_asy_airy_vpa(n_num, theta);  % expect vpa outputs
        dtheta = -val0 ./ ( vpa(sqrt(2)) * vpa(sqrt(2*n_num+1)) .* val1 .* sin(theta) );
        theta = theta - dtheta;
        if max(abs(double(dtheta))) < 1e-30  % use numeric check to avoid slow vpa norm
            break;
        end
    end

    t0 = cos(theta);
    x = vpa(sqrt(2*n_num+1) .* t0);             % back to x-variable
    w = x .* val0 + vpa(sqrt(2)) .* val1;
    w = 1 ./ (w.^2);
    %Symmetrize nodes and weights
    if mod(n,2) == 1
        w_full = [flip(w(:)); w(2:end)];
        x_full = [-flip(x(:)); x(2:end)];
    else
        w_full = [flip(w(:)); w(:)];
        x_full = [-flip(x(:)); x(:)];
    end

    w_full = vpa(w_full * sqrt(sym(pi)) / sum(exp(-vpa(x_full).^2) .* vpa(w_full)));

    x = x_full;
    w = w_full;
    % normalize weights (use vpa arithmetic)
    % if mod(n_num,2) == 1
    %     s = 2*sum(exp(-x.^2).*w) - w(1);
    % else
    %     s = 2*sum(exp(-x.^2).*w);
    % end
    % w = w * vpa(sqrt(pi)) / s;
end
  function x_init = hermite_initialguess_num(n_num)
    % numeric double initial guesses for Hermite nodes (Gatteschi + Tricomi)
    if n_num < 20
        error('n must be >= 20 for asymptotic initial guesses');
    end

    if mod(n_num,2)==1
        m = (n_num-1)/2;
        a = 0.5;
    else
        m = n_num/2;
        a = -0.5;
    end
    nu = 4*m + 2*a + 2;

    % AIRY-based roots (numeric)
    idx = (1:m)';
    T = @(t) t.^(2/3).*(1 + 5/48*t.^(-2) - 5/36*t.^(-4) + ...
        (77125/82944)*t.^(-6) - 108056875/6967296*t.^(-8) + ...
        162375596875/334430208*t.^(-10));

    airyrts = -T(3*pi/8*(4*idx-1));
    airyrts(1:min(10,m)) = vpa([-2.33810741045976703848919725244674, -4.08794944413097061663698870145739, -5.52055982809555105912985551293129, ...
                       -6.78670809007175899878024638449618, -7.94413358712085312313828055579827, -9.02265085334098038015819083988009, ...
                       -10.04017434155808593059455673736252, -11.00852430373326289323543964959015, -11.93601556323626251700636490293058, ...
                       -12.82877675286575720040672940724182]);
    x_init_airy = sqrt(abs(nu + (2^(2/3))*airyrts*nu^(1/3) + ...
        (1/5*2^(4/3))*airyrts.^2*nu^(-1/3) + ...
        (11/35-a^2-12/175)*airyrts.^3/nu));
    x_init_airy = flip(x_init_airy(:));

    % Tricomi near zero (numeric)
    Tnk0 = (pi/2)*ones(m,1);
    rhs = ((4*m+3) - 4*(1:m))/nu * pi;
    for k = 1:7
        val = Tnk0 - sin(Tnk0) - rhs(:);
        dval = 1 - cos(Tnk0);
        Tnk0 = Tnk0 - val./dval;
    end
    tnk = cos(Tnk0/2).^2;
    x_init_sin = sqrt(nu*tnk - ((tnk+1/4)./(tnk-1).^2 + (3*a^2-1))/(3*nu));
    % patch
    p = 0.4985 + eps;
    idx1 = floor(p*n_num);
    idx2 = ceil(p*n_num);
    x_init = [x_init_sin(1:idx1); x_init_airy(idx2:end)];

    if mod(n_num,2)==1
        x_init = [0; x_init];
        x_init = x_init(1:(m+1));
    else
        x_init = x_init(1:m);
    end
  end
function [val, dval] = hermpoly_asy_airy_vpa(n_num, theta)
    % theta is vpa array
    % Precompute vpa constants
    PIV = vpa(pi);
    musq = vpa(2*n_num + 1);

    cosT = cos(theta);
    sinT = sin(theta);
    sin2T = 2*cosT.*sinT;

    eta = vpa(0.5)*theta - vpa(0.25)*sin2T;
    chi = -( (vpa(3)*eta/2) ).^(vpa(2)/vpa(3));  % (-(3*eta/2)^(2/3))
    % careful with signs/powers; ensure vpa powers
    phi = (-chi ./ sinT.^2).^(vpa(1)/vpa(4));   % 1/4

    C = 2*sqrt(vpa(pi))*musq^(vpa(1)/vpa(6)).*phi;
    z = musq.^(vpa(2)/vpa(3)).*chi;

    % Convert z to double for MATLAB airy, but clip if too large
    z_num = double(z);
    Ai_num = real(airy(0, z_num));
    Aip_num = real(airy(1, z_num));
    Ai = vpa(Ai_num);
    Aip = vpa(Aip_num);

    % u-polynomials (use vpa)
    u0 = vpa(1);
    u1 = (cosT.^3 - 6*cosT)/24;
    u2 = (-9*cosT.^4 + 249*cosT.^2 + 145)/1152;
    u3 = (-4042*cosT.^9 + 18189*cosT.^7 - 28287*cosT.^5 ...
        - 151995*cosT.^3 - 259290*cosT)/414720;

    % coefficients (a,b) as vpa
    a0=vpa(1); b0=vpa(1);
    a1=vpa(15)/vpa(144); b1=-vpa(7)/vpa(5)*a1;
    a2=vpa(5*7*9*11)/vpa(2)/vpa(144)^2; b2=-vpa(13)/vpa(11)*a2;
    a3=vpa(7*9*11*13*15*17)/vpa(6)/vpa(144)^3; b3=-vpa(19)/vpa(17)*a3;

    % Hermite polynomial value (vpa arithmetic)
    val = Ai;
    B0 = -(a0*u1.*phi.^6 + a1*u0) ./ chi.^2;
    val = val + B0.*Aip/musq^(vpa(4)/vpa(3));

    A1 = (b0*u2.*phi.^12 + b1*u1.*phi.^6 + b2*u0) ./ chi.^3;
    val = val + A1.*Ai/musq^2;

    B1 = -(u3.*phi.^18 + a1*u2.*phi.^12 + a2*u1.*phi.^6 + a3*u0) ./ chi.^5;
    val = val + B1.*Aip/musq^(vpa(4)/vpa(3)+2);

    val = C.*val;

    % Derivative part
    C = sqrt(2*vpa(pi)).*musq.^(vpa(1)/vpa(3))./phi;

    v0 = vpa(1);
    v1 = (cosT.^3 + 6*cosT)/24;
    v2 = (15*cosT.^4 - 327*cosT.^2 - 143)/1152;
    v3 = (259290*cosT + 238425*cosT.^3 - 36387*cosT.^5 ...
        + 18189*cosT.^7 - 4042*cosT.^9)/414720;

    C0 = -(b0*phi.^6.*v1 + b1*v0)./chi;
    dval = C0.*Ai/musq.^(vpa(2)/vpa(3));

    D0 = a0*v0;
    dval = dval + D0.*Aip;

    C1 = -(v3.*phi.^18 + b1*v2.*phi.^12 + b2*v1.*phi.^6 + b3*v0)./chi.^4;
    dval = dval + C1.*Ai/musq.^(vpa(2)/vpa(3)+2);

    D1 = (a0*v2.*phi.^12 + a1*v1.*phi.^6 + a2*v0)./chi.^3;
    dval = dval + D1.*Aip/musq^2;

    dval = C.*dval;



end
