clear; clc;

%% =============================================================
% Pure MATLAB mesh: circular ABC + square scatterer (constrained)
% 목표: 삼각(육각) 격자 기반 + 원/사각 경계 모두 제약(CDT) + 내부만 유지
%% =============================================================
rABC = 5.0;      % outer circular ABC radius
aSq  = 2.0;      % inner square side
h    = 0.06;     % target spacing

tolUniq = 1e-12;     % uniquetol tolerance
eps0    = 1e-12;

%% -------------------------------------------------------------
% 1) Background nodes (triangular lattice) inside circle
%% -------------------------------------------------------------
dy = h*sqrt(3)/2;
ys = (-rABC:dy:rABC).';
pts = zeros(0,2);

for iy = 1:numel(ys)
    y0 = ys(iy);
    xoff = mod(iy,2) * (h/2);
    xs = (-rABC:h:rABC).' + xoff;
    row = [xs, y0*ones(size(xs))];
    pts = [pts; row];
end

mask = (pts(:,1).^2 + pts(:,2).^2) <= (rABC^2 + 10*eps0);
pts  = pts(mask,:);

%% -------------------------------------------------------------
% 2) Square boundary nodes (CCW, corners included exactly once)
%% -------------------------------------------------------------
nsq = max(4, round(aSq/h));
t   = linspace(-aSq/2, aSq/2, nsq+1).';

bottom = [ t(1:end-1), -aSq/2*ones(nsq,1) ];                % (-1,-1) 포함, ( 1,-1) 제외
right  = [  aSq/2*ones(nsq,1), t(1:end-1) ];                % ( 1,-1) 포함, ( 1, 1) 제외
top    = [ flipud(t(2:end)),  aSq/2*ones(nsq,1) ];          % ( 1, 1) 포함, (-1, 1) 제외
left   = [ -aSq/2*ones(nsq,1), flipud(t(2:end)) ];          % (-1, 1) 포함, (-1,-1) 제외

sqPts = [bottom; right; top; left];

%% -------------------------------------------------------------
% 3) Outer circle boundary nodes (constrained)
%% -------------------------------------------------------------
nCirc = max(64, round(2*pi*rABC/h));
th = linspace(0, 2*pi, nCirc+1).';
th(end) = [];
circPts = [rABC*cos(th), rABC*sin(th)];

%% -------------------------------------------------------------
% 4) Merge points with tolerance + build constraints
%% -------------------------------------------------------------
ptsAll_raw = [sqPts; circPts; pts];

[ptsAll, ~, ic] = uniquetol(ptsAll_raw, tolUniq, 'ByRows', true);

nSq   = size(sqPts,1);
nCirc2 = size(circPts,1);

idxSq   = ic(1:nSq);
idxCirc = ic(nSq + (1:nCirc2));

idxSq   = idxSq([true; diff(idxSq)~=0]);
idxCirc = idxCirc([true; diff(idxCirc)~=0]);

C_sq   = [idxSq,   idxSq([2:end 1])];
C_circ = [idxCirc, idxCirc([2:end 1])];

C = [C_sq; C_circ];
C(C(:,1)==C(:,2),:) = [];


%% -------------------------------------------------------------
% 5) Constrained Delaunay triangulation + keep triangles inside circle
%% -------------------------------------------------------------
DT = delaunayTriangulation(ptsAll, C);

P_all = DT.Points;
T_all = DT.ConnectivityList;

% keep triangles whose centroid is inside outer circle
xc = (P_all(T_all(:,1),1) + P_all(T_all(:,2),1) + P_all(T_all(:,3),1))/3;
yc = (P_all(T_all(:,1),2) + P_all(T_all(:,2),2) + P_all(T_all(:,3),2))/3;

maskOuter = (xc.^2 + yc.^2) <= (rABC^2 + 10*eps0);

T = T_all(maskOuter,:);
P = P_all;


%% -------------------------------------------------------------
% 6) FEM format variables
%% -------------------------------------------------------------
p = P;           % nNodes x 2
N = T;           % nElem x 3

x = p(:,1);
y = p(:,2);

nNodes = size(p,1);
nElem  = size(N,1);

%% -------------------------------------------------------------
% 7) CCW orientation for each triangle
%% -------------------------------------------------------------
for e = 1:nElem
    x21 = x(N(e,2)) - x(N(e,1));
    x31 = x(N(e,3)) - x(N(e,1));
    y21 = y(N(e,2)) - y(N(e,1));
    y31 = y(N(e,3)) - y(N(e,1));
    Ae  = 0.5*(x21*y31 - x31*y21);
    if Ae < 0
        tmp   = N(e,2);
        N(e,2)= N(e,3);
        N(e,3)= tmp;
    end
end

%% -------------------------------------------------------------
% 8) ABC boundary extraction (use constrained circle edges)
%% -------------------------------------------------------------
TR = triangulation(N, p);
F  = freeBoundary(TR);

rnode = hypot(p(:,1), p(:,2));
tolABC = max(2*h, 1e-3);

isABCnode   = abs(rnode - rABC) <= tolABC;
maskABCedge = isABCnode(F(:,1)) & isABCnode(F(:,2));

NS      = F(maskABCedge,:);     % ABC boundary segments
ABC     = unique(NS(:));
nSegABC = size(NS,1);

%% -------------------------------------------------------------
% 9) Region labeling by centroid (1=background, 2=inside square)
%% -------------------------------------------------------------
xc = (x(N(:,1)) + x(N(:,2)) + x(N(:,3))) / 3;
yc = (y(N(:,1)) + y(N(:,2)) + y(N(:,3))) / 3;

region = ones(nElem,1);
tolReg = 1e-9;
insideElem = (abs(xc) < aSq/2 - tolReg) & (abs(yc) < aSq/2 - tolReg);
region(insideElem) = 2;
%% =============================================================
% Build ElemOfSeg for boundary segments NS (nSegABC x 2)
% Input: N (nElem x 3), NS (nSegABC x 2), nNodes
% Output: ElemOfSeg (nSegABC x 1)
%% =============================================================

nNodes   = size(p,1);
nElem    = size(N,1);
nSegABC  = size(NS,1);

% 1) 모든 요소 edge 3개씩 뽑기 (총 3*nElem)
E12 = N(:,[1 2]);
E23 = N(:,[2 3]);
E31 = N(:,[3 1]);

Eall = [E12; E23; E31];                % (3*nElem) x 2
Eall = sort(Eall, 2);                  % [min max] 로 정규화

elem_id = repmat((1:nElem).', 3, 1);   % 각 edge가 속한 요소 id

% 2) edge를 유일키로 변환: key = i + (j-1)*nNodes  (i<j)
key_all = Eall(:,1) + (Eall(:,2)-1)*nNodes;

% 경계 edge는 요소 1개만 가져야 정상.
% 혹시 같은 edge가 여러 번 들어오면(비정상), 첫 번째만 사용.
[key_u, ia] = unique(key_all, 'stable');
elem_u      = elem_id(ia);

% 3) NS도 같은 방식으로 key 만들고 매칭
NSs   = sort(NS, 2);
key_NS = NSs(:,1) + (NSs(:,2)-1)*nNodes;

[tf, loc] = ismember(key_NS, key_u);

if any(~tf)
    bad = find(~tf);
    error('ElemOfSeg 매칭 실패: NS의 %d개 edge가 요소 edge 목록에 없습니다. 예: s=%d (nodes %d-%d)', ...
          nnz(~tf), bad(1), NS(bad(1),1), NS(bad(1),2));
end

ElemOfSeg = elem_u(loc);   % (nSegABC x 1)
%% -------------------------------------------------------------
% 10) Debug plot
%% -------------------------------------------------------------
figure; hold on; axis equal;
triplot(N, x, y, 'k');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 1;
E_0 = 1; % Amplitude of incident electric field (V/m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE CONSTANTS FOR EACH ELEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J = complex(0.,1.);
eps_sq = 3.0 - 0.5*J;    % 정사각형 내부
eps_bg = 1.0;           % 도넛 영역
mu_r   = 1.0;
k_0 = 2*pi/lambda; % wavenumber
alpha = -1/mu_r;
kappa = 1/rABC;
gamma_1 = alpha * ( J*k_0 + kappa/2 ...
                    - J*kappa^2/(8*(J*kappa - k_0)) );

gamma_2 = -J*alpha / (2*(J*kappa - k_0));
theta = 0; 
mu_r = zeros(1,nElem); % relative permeability
eps_r = zeros(1,nElem); % relative permittivity
alpha_x = zeros(1,nElem); % coefficient of second order x term in PDE
alpha_y = zeros(1,nElem); % coefficient of second order y term in PDE
beta = zeros(1,nElem); % coefficient of term linear in field in PDE
A = zeros(1,nElem); % area of element e

for e = 1:nElem
    if region(e) == 1
        eps_r(e) = eps_bg;
        mu_r(e)  = 1.0;
    else
        eps_r(e) = eps_sq;
        mu_r(e)  = 1.0;
    end

    alpha_x(e) = -1/mu_r(e);
    alpha_y(e) = -1/mu_r(e);
    beta(e)    = (k_0^2) * eps_r(e);   % 복소수
end

% Initialize global K matrix and the right hand side vector b
Kglobal = zeros(nNodes, nNodes);
btilde = zeros(nNodes, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMATION OF ELEMENT MATRICES AND ASSEMBLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 선형 삼각형 보간함수 계수
a = zeros(nElem,3);
b = zeros(nElem,3);
c = zeros(nElem,3);

be_tilde = zeros(nElem,3);   % 요소별 RHS
Klocal   = zeros(nElem,3,3); % 요소 stiffness+mass

for e = 1:nElem
    % --- 좌표 및 면적 ---
    x21 = x(N(e,2)) - x(N(e,1));
    x31 = x(N(e,3)) - x(N(e,1));
    y21 = y(N(e,2)) - y(N(e,1));
    y31 = y(N(e,3)) - y(N(e,1));

    A(e) = 0.5*(x21*y31 - x31*y21);   % 삼각형 면적
    Ae = abs(A(e));
    
    x1 = x(N(e,1)); y1 = y(N(e,1));
    x2 = x(N(e,2)); y2 = y(N(e,2));
    x3 = x(N(e,3)); y3 = y(N(e,3));

    % --- shape function 계수 a,b,c ---
    a(e,1) = x2*y3 - x3*y2;
    b(e,1) = y2 - y3;
    c(e,1) = x3 - x2;

    a(e,2) = x3*y1 - x1*y3;
    b(e,2) = y3 - y1;
    c(e,2) = x1 - x3;

    a(e,3) = x1*y2 - x2*y1;
    b(e,3) = y1 - y2;
    c(e,3) = x2 - x1;

    % --- 요소 K^e 행렬 ---
    for i = 1:3
        for j = 1:3
            if i == j
                delta_ij = 1;
            else
                delta_ij = 0;
            end

            Klocal(e,i,j) = ( ...
                alpha_x(e)*b(e,i)*b(e,j) + ...
                alpha_y(e)*c(e,i)*c(e,j) )/(4*Ae) ...
                + Ae*beta(e)*(1 + delta_ij)/12;
        end
    end

    % --- 요소 b_tilde^e 벡터 (산란장 공식) ---
    for i = 1:3
        fun = @(u,v) betilde(u,v, ...
            a(e,i), b(e,i), c(e,i), ...
            k_0, x1, x2, x3, y1, y2, y3,theta);

        be_tilde(e,i) = (E_0*k_0*k_0*(1/mu_r(e) - eps_r(e))/(2*Ae)) * ...
                        dblquad(fun,0,1,0,1);   % 필요하면 integral2로 교체 가능
    end
end

% ---------------- 전역 행렬/벡터 조립 ----------------
for e = 1:nElem
    for i = 1:3
        Ii = N(e,i);
        for j = 1:3
            Jj = N(e,j);
            Kglobal(Ii,Jj) = Kglobal(Ii,Jj) + Klocal(e,i,j);
        end
        btilde(Ii) = btilde(Ii) + be_tilde(e,i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ABC 경계(2차) 행렬/벡터 (Ks, bs_tilde) 생성 및 조립
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ks       = zeros(nSegABC,2,2);
bs_tilde = zeros(nSegABC,2);

for s = 1:nSegABC
    n1 = NS(s,1);  n2 = NS(s,2);
    x1 = x(n1);    x2 = x(n2);
    y1 = y(n1);    y2 = y(n2);

    L_s = hypot(x1-x2, y1-y2);

    for ii = 1:2
        % ① 1차 → 2차: bstilde1 대신 bstilde2 호출
        fun2 = @(ksi) bstilde2(ksi, ii, alpha, kappa, k_0, ...
                               x1, x2, y1, y2, L_s, E_0, theta);
        bs_tilde(s,ii) = quad(fun2, 0, 1);

        % ② Ks 구성식 변경: gamma_2 / L_s 항 추가
        for jj = 1:2
            if ii == jj
                Ks(s,ii,jj) = gamma_1*L_s/3 - gamma_2 / L_s;
            else
                Ks(s,ii,jj) = gamma_1*L_s/6 + gamma_2 / L_s;
            end
        end
    end
end


% --- ABC 경계 행렬/벡터를 전역에 조립 ---
for s = 1:nSegABC
    for i = 1:2
        Ii = NS(s,i);
        for j = 1:2
            Jj = NS(s,j);
            Kglobal(Ii,Jj) = Kglobal(Ii,Jj) + Ks(s,i,j);
        end
        btilde(Ii) = btilde(Ii) + bs_tilde(s,i);
    end
end

%% -------------------------------------------------------------
%  Solution of global matrix system
%% -------------------------------------------------------------
Ez = Kglobal \ btilde;   % 복소 스칼라장 Ez (nNodes x 1)

%% -------------------------------------------------------------
%  Generate solution over a grid and interpolate FEM solution
%% -------------------------------------------------------------
[xgrid, ygrid] = meshgrid(-rABC : 0.01*(2*rABC) : rABC, ...
                          -rABC : 0.01*(2*rABC) : rABC);

[nx, ny] = size(xgrid);
Ezgrid   = zeros(nx, ny);    % 격자 위의 Ez (복소 또는 실수)

for i = 1:nx
    for j = 1:ny

        xv = xgrid(i,j);
        yv = ygrid(i,j);

        % 이 (xv,yv)가 포함된 요소 e를 찾고, 그 요소의 shape function으로 보간
        for e = 1:nElem

            % --- 점 (xv,yv)가 삼각형 e 안에 있는지 검사하기 위한 부분면적법 ---
            x2p = x(N(e,2)) - xv;
            x3p = x(N(e,3)) - xv;
            y2p = y(N(e,2)) - yv;
            y3p = y(N(e,3)) - yv;
            A1  = 0.5*abs(x2p*y3p - x3p*y2p);

            x2p = x(N(e,2)) - xv;
            x1p = x(N(e,1)) - xv;
            y2p = y(N(e,2)) - yv;
            y1p = y(N(e,1)) - yv;
            A2  = 0.5*abs(x2p*y1p - x1p*y2p);

            x1p = x(N(e,1)) - xv;
            x3p = x(N(e,3)) - xv;
            y1p = y(N(e,1)) - yv;
            y3p = y(N(e,3)) - yv;
            A3  = 0.5*abs(x1p*y3p - x3p*y1p);

            x21 = x(N(e,2)) - x(N(e,1));
            x31 = x(N(e,3)) - x(N(e,1));
            y21 = y(N(e,2)) - y(N(e,1));
            y31 = y(N(e,3)) - y(N(e,1));

            % Ae는 부호 있는 면적(원래 FEM 형식 그대로 사용하는 게 안전)
            Ae  = 0.5*(x21*y31 - x31*y21);

            % 세 부분면적 합이 전체면적과 거의 같으면, (xv,yv)는 삼각형 내부
            if abs(Ae - (A1 + A2 + A3)) < 1.0e-5*abs(Ae)

                % --- barycentric 좌표(ksi, ita) 계산 ---
                ksi = ( y31*(xv - x(N(e,1))) - x31*(yv - y(N(e,1))) ) / (2*Ae);
                ita = (-y21*(xv - x(N(e,1))) + x21*(yv - y(N(e,1))) ) / (2*Ae);

                % 선형 shape function
                N1 = 1 - ksi - ita;
                N2 = ksi;
                N3 = ita;

                % FEM 노드 값으로부터 Ez(xv,yv) 보간
                Ezgrid(i,j) = N1*Ez(N(e,1)) + N2*Ez(N(e,2)) + N3*Ez(N(e,3));

                break;   % 해당 요소 찾았으니 e-loop 탈출
            end
        end
    end
end

%% -------------------------------------------------------------
%  Display contour plot of FEM solution
%% -------------------------------------------------------------
figure;
contourf(xgrid, ygrid, real(Ezgrid), 20, 'LineColor', 'none');
xlabel('x (wavelengths)');
ylabel('y (wavelengths)');
axis([-rABC rABC -rABC rABC]);
axis square;
title('Total Electric Field |E_z| - Contour Plot');
colorbar;

%% -------------------------------------------------------------
%  Evaluate FEM field Ez on a circle r = dist (near-field cut)
%% -------------------------------------------------------------
d2p = pi/180;

dist = 2.5*lambda;                 % 원하는 평가 반경 (예: lambda, 2*lambda 등)
dphi_deg = 0.25;               % 각도 샘플링 간격 (당신 코드: 0.25 deg)
phi_deg = 0:dphi_deg:360;      % 0~360도
nPhi = numel(phi_deg);

Ez_eval = nan(1,nPhi);

for I = 1:nPhi
    ph = phi_deg(I)*d2p;
    xeval = dist*cos(ph);
    yeval = dist*sin(ph);

    % 점이 속한 삼각형을 찾아 보간
    for e = 1:nElem
        n1 = N(e,1); n2 = N(e,2); n3 = N(e,3);

        x1 = x(n1); y1 = y(n1);
        x2 = x(n2); y2 = y(n2);
        x3 = x(n3); y3 = y(n3);

        % 부호 있는 면적 (orientation 따라 +/-)
        Ae_signed = 0.5*((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
        if Ae_signed == 0
            continue;
        end
        Ae = abs(Ae_signed);

        % 부분면적(절댓값)으로 내부판정
        A1 = 0.5*abs( (x2-xeval)*(y3-yeval) - (x3-xeval)*(y2-yeval) );
        A2 = 0.5*abs( (x3-xeval)*(y1-yeval) - (x1-xeval)*(y3-yeval) );
        A3 = 0.5*abs( (x1-xeval)*(y2-yeval) - (x2-xeval)*(y1-yeval) );

        if abs(Ae - (A1 + A2 + A3)) < 1.0e-6*Ae
            % barycentric 좌표 (orientation-safe하게 계산)
            % l1,l2,l3 such that r = l1*r1 + l2*r2 + l3*r3, l1+l2+l3=1
            denom = ( (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3) );  % = 2*Ae_signed
            if denom == 0
                break;
            end
            l1 = ( (y2-y3)*(xeval-x3) + (x3-x2)*(yeval-y3) ) / denom;
            l2 = ( (y3-y1)*(xeval-x3) + (x1-x3)*(yeval-y3) ) / denom;
            l3 = 1 - l1 - l2;

            Ez_eval(I) = l1*Ez(n1) + l2*Ez(n2) + l3*Ez(n3);
            break;
        end
    end
end

% 혹시 원이 도메인 밖이면 nan이 생깁니다.
if any(isnan(Ez_eval))
    warning('Ez_eval has NaNs. Check if dist is inside the FEM computational domain (inside ABC circle).');
end

% Plot
figure;
plot(phi_deg, abs(Ez_eval), 'b-');
grid on;
xlabel('\phi (deg)');
ylabel('|E_z(\rho=dist,\phi)|');
title(sprintf('FEM field magnitude on circle r = %.3g', dist));

% %% ============================================
% %  Field sampling on a circle: Ez(phi) at r=dist
% %  (FEM nodal solution Ez -> interpolate on circle)
% %% ============================================
% d2p  = pi/180;
% k0   = k_0;
% 
% % NTFF / RCS evaluation radius (rho in Jin)
% dist = 2.0*lambda;
% 
% % incidence angles
% theta_deg = (0:45).';
% theta_rad = theta_deg*d2p;
% nTheta    = numel(theta_rad);
% 
% nNodes = numel(x);
% nElem  = size(N,1);
% nSegABC = size(NS,1);
% 
% 
% 
% % --- factorize once ---
% [L,U,Pperm] = lu(Kglobal);
% 
% % --- solve u for all incidence angles ---
% u = zeros(nNodes, nTheta);
% for it = 1:nTheta
%     theta = theta_rad(it);
% 
%     % ==============================
%     % (NEW) Volume RHS: angle-dependent
%     % ==============================
%     btilde_vol = zeros(nNodes,1);
% 
%     for e = 1:nElem
%         n1e = N(e,1); n2e = N(e,2); n3e = N(e,3);
% 
%         x1e = x(n1e); y1e = y(n1e);
%         x2e = x(n2e); y2e = y(n2e);
%         x3e = x(n3e); y3e = y(n3e);
%         pref = (E_0*k0^2*(1/mu_r(e) - eps_r(e)));
% 
%         for i = 1:3
%             funV = @(u,v) betilde_tri(u,v, a(e,i), b(e,i), c(e,i), ...
%                          k0, x1e, y1e, x2e, y2e, x3e, y3e, theta);
% 
%             be_i = pref * integral2(funV, ...
%             0, 1, ...
%             @(u) zeros(size(u)), ...   % 중요: u와 같은 크기
%             @(u) 1 - u, ...
%             'AbsTol',1e-9,'RelTol',1e-6);
% 
% 
%             Ii = N(e,i);
%             btilde_vol(Ii) = btilde_vol(Ii) + be_i;
%         end
%     end
% 
%     % ==============================
%     % ABC RHS (angle-dependent) - keep as you have
%     % ==============================
%     btilde_ang = zeros(nNodes,1);
% 
%     for s = 1:nSegABC
%         n1 = NS(s,1);
%         n2 = NS(s,2);
% 
%         x1 = x(n1);  y1 = y(n1);
%         x2 = x(n2);  y2 = y(n2);
% 
%         Ls = hypot(x2-x1, y2-y1);
% 
%         for m = 1:2
%             funB = @(ksi) bstilde1_angle(ksi, m, alpha, kappa, k0, ...
%                                         x1, y1, x2, y2, Ls, E_0, theta);
%             bs_val = quadgk(funB, 0, 1);
% 
%             Ii = NS(s,m);
%             btilde_ang(Ii) = btilde_ang(Ii) + bs_val;
%         end
%     end
% 
%     % total RHS for this theta
%     btilde = btilde_vol + btilde_ang;
% 
%     % solve
%     u(:,it) = U \ (L \ (Pperm * btilde));
% end
% 
% % --- build ElemOfSeg for boundary edges ---
% TRg = triangulation(N, [x(:) y(:)]);
% att = edgeAttachments(TRg, NS);
% ElemOfSeg = zeros(nSegABC,1);
% for s = 1:nSegABC
%     eIDs = att{s};
%     if isempty(eIDs)
%         error('No adjacent element for NS edge %d', s);
%     end
%     ElemOfSeg(s) = eIDs(1);
% end
% 
% % --- compute RCS/λ in dB via Jin (1.88)-(1.89) ---
% RCS_over_lambda_dB = zeros(nTheta,1);
% 
% for it = 1:nTheta
%     theta = theta_rad(it);
% 
%     Esc = NTFF_Esc_Jin(u(:,it), theta, x, y, NS, ElemOfSeg, N, A, k0, E_0, dist);
% 
%     sigma = 2*pi*dist * (abs(Esc)^2) / (abs(E_0)^2);     % Jin (1.88)
%     RCS_over_lambda_dB(it) = 10*log10(sigma/lambda);     % σ/λ in dB
% end
% 
% % --- plot ---
% figure;
% plot(theta_deg, RCS_over_lambda_dB, 'b-'); grid on;
% xlabel('\theta^{inc} (deg)');
% ylabel('10log_{10}(\sigma/\lambda) (dB)');
% title('Monostatic RCS/\lambda (Jin 1.88-1.89)');
% 
% %% =========================
% %% FUNCTIONS (paste below)
% %% =========================
% function f1 = bstilde1_angle(ksi, m, alpha, kappa, k0, ...
%                              x1, y1, x2, y2, Ls, E0, theta)
%     if m == 1
%         Nsh = 1 - ksi;
%     else
%         Nsh = ksi;
%     end
% 
%     x = (1-ksi).*x1 + ksi.*x2;
%     y = (1-ksi).*y1 + ksi.*y2;
% 
%     phase = x*cos(theta) + y*sin(theta);
%     expfac = exp(-1i*k0.*phase);
% 
%     coeff = (1i*k0 + 0.5*kappa);
%     f1 = alpha .* coeff .* Nsh .* Ls .* E0 .* expfac;
% end
% 
% function Esc = NTFF_Esc_Jin(u, theta, x, y, NS, ElemOfSeg, N, A, k0, E0, rho)
%     % monostatic observation direction
%     rhat = [-cos(theta); -sin(theta)];
% 
%     % 2-pt Gauss on [0,1]
%     ksi = [0.211324865405187, 0.788675134594813];
%     w   = [0.5, 0.5];
% 
%     I = 0;
% 
%     for s = 1:size(NS,1)
%         n1 = NS(s,1); n2 = NS(s,2);
% 
%         x1 = x(n1); y1 = y(n1);
%         x2 = x(n2); y2 = y(n2);
% 
%         Ls = hypot(x2-x1, y2-y1);
% 
%         % outward normal for circular ABC (orientation-independent)
%         xm = 0.5*(x1+x2); ym = 0.5*(y1+y2);
%         rm = hypot(xm,ym);
%         nx = xm/rm; ny = ym/rm;
% 
%         e = ElemOfSeg(s);
%         grad_utot = grad_u_element_P1(e, u, x, y, N, A);
% 
%         for q = 1:2
%             xp = (1-ksi(q))*x1 + ksi(q)*x2;
%             yp = (1-ksi(q))*y1 + ksi(q)*y2;
% 
%             u_tot = (1-ksi(q))*u(n1) + ksi(q)*u(n2);
% 
%             [u_inc, duninc_dn] = incident_plane_wave(xp, yp, nx, ny, k0, E0, theta);
% 
%             u_s = u_tot - u_inc;
%             dun_s_dn = (grad_utot(1)*nx + grad_utot(2)*ny) - duninc_dn;
% 
%             phase = exp(1i*k0*(rhat(1)*xp + rhat(2)*yp));
% 
%             % Jin (1.89) integrand: [ jk0(rhat·n) u_s - ∂u_s/∂n ] e^{jk0 rhat·r'}
%             I = I + ( 1i*k0*(rhat(1)*nx + rhat(2)*ny)*u_s - dun_s_dn ) * phase * w(q) * Ls;
%         end
%     end
% 
%     % Jin (1.89) prefactor
%     Esc = sqrt(1i*k0/(8*pi*rho)) * exp(-1i*k0*rho) * I;
% end
% 
% function grad_u = grad_u_element_P1(e, u, x, y, N, A)
%     n1 = N(e,1); n2 = N(e,2); n3 = N(e,3);
% 
%     x1 = x(n1); y1 = y(n1);
%     x2 = x(n2); y2 = y(n2);
%     x3 = x(n3); y3 = y(n3);
% 
%     u1 = u(n1); u2 = u(n2); u3 = u(n3);
%     Ae = A(e);
% 
%     gradN1 = [ y2 - y3 ; x3 - x2 ] / (2*Ae);
%     gradN2 = [ y3 - y1 ; x1 - x3 ] / (2*Ae);
%     gradN3 = [ y1 - y2 ; x2 - x1 ] / (2*Ae);
% 
%     grad_u = u1*gradN1 + u2*gradN2 + u3*gradN3;
% end
% 
% function [uinc, duninc_dn] = incident_plane_wave(xp, yp, nx, ny, k0, E0, theta)
%     phase = xp*cos(theta) + yp*sin(theta);
%     uinc = E0 * exp(-1i*k0*phase);
% 
%     duninc_dn = (-1i*k0*uinc) * (cos(theta)*nx + sin(theta)*ny);
% end
% 
% function val = betilde_tri(u,v, ai, bi, ci, k0, x1,y1, x2,y2, x3,y3, theta)
%     % Reference triangle: (u,v), u>=0, v>=0, u+v<=1
%     % Physical mapping: r = r1 + u*(r2-r1) + v*(r3-r1)
% 
%     x = x1 + u.*(x2-x1) + v.*(x3-x1);
%     y = y1 + u.*(y2-y1) + v.*(y3-y1);
% 
%     % Constant Jacobian magnitude for affine triangle mapping
%     J = abs( (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) );  % = 2*Area
% 
%     phase = x*cos(theta) + y*sin(theta);
% 
%     % P1 basis function for node i: Ni = (ai + bi*x + ci*y) / (2A)
%     % but since J = 2A, we can write Ni * dA = (ai + bi*x + ci*y)/J * (J du dv)
%     % => (ai + bi*x + ci*y) * du dv
%     % Therefore integrand is simply (ai + bi*x + ci*y) * exp(-j k0 phase).
%     val = (ai + bi.*x + ci.*y) .* exp(-1i*k0.*phase);
% end
