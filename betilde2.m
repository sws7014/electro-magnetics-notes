function f1 = bstilde2(ksi,m,alpha,kappa,k_0,...
                       x1,x2,y1,y2,L_s,E_0,theta)

    if m == 1
        N = 1-ksi;
    else
        N = ksi;
    end

    J = 1i;

    x = (1-ksi).*x1 + ksi.*x2;
    y = (1-ksi).*y1 + ksi.*y2;

    % rotated coordinates
    s = x*cos(theta) + y*sin(theta);        % propagation
    t = -x*sin(theta) + y*cos(theta);       % transverse

    expfac = exp(-J*k_0.*s);

    term1 = ( J*k_0 + 0.5*kappa ...
              - J*kappa^2/(8*(J*kappa - k_0)) );

    term2 = ( - J*k_0*kappa.*s ...
              - 0.5*J*kappa^2/(J*kappa - k_0) .* ...
                ( J*k_0.*s - k_0^2.*t.^2 ) );

    f1 = alpha * ( term1 + term2 ) .* N .* L_s .* E_0 .* expfac;
end
