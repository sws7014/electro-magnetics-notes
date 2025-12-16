function f2 = betilde(u,v,a,b,c,k_0,x1,x2,x3,y1,y2,y3,theta)
J = 1i;

x = (1-u).*x1 + u.*((1-v).*x2 + v.*x3);
y = (1-u).*y1 + u.*((1-v).*y2 + v.*y3);

jacob = (-x1+(1-v).*x2+v.*x3).*(-u.*y2+u.*y3) ...
      - (-y1+(1-v).*y2+v.*y3).*(-u.*x2+u.*x3);

phase = x*cos(theta) + y*sin(theta);

f2 = (a+b.*x+c.*y).*exp(-J*k_0.*phase).*jacob;
end
