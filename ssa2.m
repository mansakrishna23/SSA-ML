% 
% EARS 107 Final Project - Mansa Krishna
%
% Shallow Shelf Approximation (2D)
% Finite Element Method
%
% Code was adapted from the following resources: 
% Larson, M. G., & Bengzon, F. (2013). The Finite Element Method: Theory, 
%   Implementation, and Applications (Vol. 10). Berlin, Heidelberg: Springer 
%   Berlin Heidelberg. https://doi.org/10.1007/978-3-642-33287-6
% Mathieu Morlighem & Helene Seroussi. (2018, October 12). Finite Element 
%   Formumation of the Shelfy Stream Approximation (SSA).

% Physical Parameters (kept constant for simplicity)
H   = 500; % next steps: set up a linear height profile
mu  = 1000;
g   = 9.81;
rho = 917;

% Minimum x and y positions
min = [0, 0];
% Maximum x and y positions
max = [10, 10];
% Set the domain using helper function
domain = SetDomain(min, max); 
% point matrix, element matrix, connectivity matrix
maxh = 1; % Maximum height of each element, maxh 
[p, e, t] = initmesh(domain, 'hmax', maxh);

tiledlayout(2,2);
nexttile
% pdeplot(p, e, t, "NodeLabels", 'on', 'ElementLabels', 'on');
pdeplot(p, e, t);
title("Meshed Domain"); xlabel("x (m)"); ylabel("y (m)");
% saveas(gcf,'mesheg.png') % save figure

% Getting the number of nodes and elements
nods = size(p, 2); 
nels = size(t, 2);

% Global stiffness matrix, K
K = zeros(2*nods, 2*nods);

for T = 1:nels
    % loop through each of the elements
    idx = t(1:3, T); % element to global index
    % x and y positions of the nodes of current element
    x = p(1, idx);
    y = p(2, idx); 
    % area of element
    area = polyarea(x, y);
    A2 = area^2; % square of the area of the element

    % Shape function (phi)i = ai + bix + ciy
    % Computing coefficients b and c
    b = [(y(2)-y(3))/(2*area)
        (y(3)-y(1))/(2*area)
        (y(1)-y(2))/(2*area)];
    c = [(x(3)-x(2))/(2*area)
        (x(1)-x(3))/(2*area)
        (x(2)-x(1))/(2*area)];

    Ke = zeros(6, 6);
    % Assembling the element stiffness matrix, Ke
    for i=1:3
        for j=1:3
            Ke(2*i-1, 2*j-1) = A2*(4*H*mu*b(j)*b(i) + H*mu*c(j)*c(i));
            Ke(2*i-1, 2*j) = A2*(2*H*mu*c(j)*b(i) + H*mu*b(j)*c(i));
            Ke(2*i, 2*j-1) = A2*(2*H*mu*b(j)*c(i) + H*mu*c(j)*b(i));
            Ke(2*i, 2*j) = A2*(4*H*mu*c(j)*c(i) + H*mu*b(j)*b(i));
        end
    end

    % indices for the global stiffness matrix
    gidx = [2*idx(1)-1; 2*idx(1); 2*idx(2)-1; 2*idx(2); 2*idx(3)-1; 2*idx(3)];
    % adding element stiffness matrix to global matrix
    K(gidx, gidx) = K(gidx, gidx) + Ke;

end

% Assembling the load vector, F
F = zeros(2*nods, 1);
for T = 1:nels
    idx = t(1:3, T); % element to global 
    x = p(1, idx); y = p(2, idx); % x and y values of nodes

    % Area of element
    area = polyarea(x,y);
    A2 = area^2; % square of the area of the element

    % surface slopes
    % try playing around with these slopes!
    dsdx = 0;
    dsdy = -0.0001;

    % Element load vector, Fe
    I = 1:1:3; % [1, 2, 3]
    Fe = zeros(6, 1);
    % Assembling Fe
    Fe(2*I-1) = -rho*g*H*(area/3)*A2*dsdx;
    Fe(2*I) = -rho*g*H*(area/3)*A2*dsdy;

    % global index
    gidx = [2*idx(1)-1; 2*idx(1); 2*idx(2)-1; 2*idx(2); 2*idx(3)-1; 2*idx(3)];
    % Adding the element load vector to the global load vector
    F(gidx) = F(gidx) + Fe;

end

% Enforcing boundary conditions
boundary = unique([e(1,:), e(2,:)]);

for i=1:numel(boundary)
    % degrees of freedom in x and y
    dofx = 2*boundary(i)-1;
    dofy = 2*boundary(i);

    % Defining the boundary values within the global 
    % stiffness matrix and global load vector
    K(dofx, :) = 0;
    K(dofx, dofx) = 1e3;
    K(dofy, :) = 0;
    K(dofy, dofy) = 1e3;
    F(dofx) = 0;
    F(dofy) = 0;
end

U = K\F; % solving linear system of equations. 

% extract u and v
uidx = 1:2:(2*nods-1); % odd indices
vidx = 2:2:(2*nods); % even indices
u = U(uidx);
v = U(vidx);

% pdeplot(p,e,t, "XYData", u, "ColorMap", 'jet');
% title("v_x(x,y)")

nexttile
pdeplot(p,e,t, "flowdata", [u, v], "ColorMap", 'jet', 'Mesh', 'on');
title("Ice Flow");

% clf;

% tiledlayout(1, 2);

nexttile
pdeplot(p,e,t, "XYData", u, "ColorMap", 'jet', 'Mesh', 'on');
title("v_x(x,y)")

nexttile
pdeplot(p,e,t, "XYData", v, "ColorMap", 'jet', 'Mesh', 'on');
title("v_y(x,y)")

% saveas(gcf, 'velocityfield.png');

% --------------------------
% Helper function for Domain
% Adapted from 
% Larson, M. G., & Bengzon, F. (2013). The Finite Element Method: Theory, 
%   Implementation, and Applications (Vol. 10). Berlin, Heidelberg: Springer 
%   Berlin Heidelberg. https://doi.org/10.1007/978-3-642-33287-6
% --------------------------

% Function takes in min = [x, y] and max = [x, y] coordinates
% and sets a rectangular domain accordingly
function r = SetDomain(min, max)
r = [2 min(1) max(1) min(2) min(2) 1 0;
   2 max(1) max(1) min(2) max(2) 1 0;
   2 max(1) min(1) max(2) max(2) 1 0;
   2 min(1) min(1) max(2) min(2) 1 0]';
end
