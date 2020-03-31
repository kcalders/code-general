function [Vtot,DiamSummary] = triangulated_cylinder_surface(P,NA,H,CL,maxL,Axe,Sta)


% Triangulates the given point cloud according to a cylindrical support. 
% The cylinder support is defined by the number of sectors and the height 
% of layers. Additional info for the support is the length of each cylinder
% the maximum total length of the support. 

% Iputs
% P         Point cloud to be approximated with triangulated surface
% na        Number of sectors used in triangulation
% H         Height (in meters) of the layers used in triangulation
% cl        Cylinder length (in meters), optional, default = 1 meter
% maxL      Maximum length to be reconstructed, if not given all the point
%               cloud will be reconstructed
% Axe       Cylinder axes, optional
% Sta       Starting points of cylinders, optional

% Outputs
% Vtot      Total volume enclosed by the polyhedral surface in liters
% Diam      Diameter data of each "slice", 
%               [minimum_diam  mean_diam  maximum_diam  height_from_base] 
%               in centimeters
% Vert      Coordinates of the vertices of the polyhedron
% Facets    Index of vertices forming the triangles
% fvd       Color info of faces for patch plotting
% Normals   Normals of the triangles/facets
% Centers   Centers (coordinates) of the triangles
% Areas     Areas of the triangles

if nargin == 3
    CL = 1;
elseif nargin == 4
    maxL = 1000;
end
np = size(P,1);


%% Determine cylinders used for the support
% Divide the trunk into "cl" meter height "cylinders"

% Calculate the height of each point
if nargin == 7
    % Use the first cylinder to give the base and direction
    [~,~,h] = distances_to_line(P,Axe(1,:),Sta(1,:));
else
    % Determine the base and direction from the point cloud
    X = cov(P);
    [U,S,V] = svd(X);
    U = U(:,1);
    if U(3) < 0
        U = -U;
    end
    [~,~,h] = distances_to_line(P,U,mean(P));
end
h = h-min(h);  % zero height for the lowest point

K = max(h);
if K > maxL
    K = maxL;
end

n = ceil(K/CL);  % number of cylinders
CL = K/n;

Sta = zeros(n+1,3);
Ind = false(np,n);
for i = 1:n
    if i < n
        I = h < i*CL;
        J = h >= (i-1)*CL;
        I = I&J;
    else
        I = h < maxL;
        J = h >= (i-1)*CL;
        I = I&J;
    end
    Sta(i,:) = mean(P(I,:));
    Ind(:,i) = I;
end

I = h >= K-H/2;
Sta(n+1,:) = mean(P(I,:));
Axe = Sta(2:end,:)-Sta(1:n,:);
Len = sqrt(sum(Axe.*Axe,2));
Axe = [Axe(:,1)./Len Axe(:,2)./Len Axe(:,3)./Len];
Sta = Sta(1:n,:);
Rad = Len;

for i = 1:n
    I = Ind(:,i);
    R0 = mean(distances_to_line(P(I,:),Axe(i,:),Sta(i,:)));
    if nnz(I) > 10
        [AP, CA, R] = lscylinder(P(I,:),Sta(i,:)',Axe(i,:)',R0,0.1,0.1);
    else
        AP = Sta(i,:)';
        CA = Axe(i,:)';
        R = R0;
    end
    h = P(I,:)*CA;
    hmin = min(h);
    L = abs(max(h)-hmin);
    hpoint = CA'*AP;
    Len(i) = L;
    Sta(i,:) = AP'-(hpoint-hmin)*CA';
    Axe(i,:) = CA/norm(CA);
    Rad(i) = R;
end

if Rad(1) > 0.03
    Sta(1,:) = Sta(1,:)+H/2*Axe(1,:);
    Len(1) = Len(1)-H/2;
end

% Modify the cylinder data to form a continuous support axis
for i = 2:n
    if i < n
        Sta(i,:) = Sta(i-1,:)+Len(i-1)*Axe(i-1,:);
        V = Sta(i+1,:)-Sta(i,:);
        Axe(i,:) = V/norm(V);
        Len(i) = norm(V);
    else
        E = Sta(i,:)+Len(i)*Axe(i,:);
        Sta(i,:) = Sta(i-1,:)+Len(i-1)*Axe(i-1,:);
        V = E-Sta(i,:);
        Axe(i,:) = V/norm(V);
        Len(i) = norm(V);
    end
end


%% Determine the vertices by fitting
% Generate the refence vectors for determination of the angular sectors
U = Axe(1,:)';
[V,W] = orthonormal_vectors(U);

% Check the handedness of the reference
if det([U V W]) < 0
    v = V;
    V = W;
    W = v;
end

% Define the vertices of the polyhedra
S = 1:1:np;
a = ceil(n*NA*CL/H)+10*n;
Vert = zeros(a,3);  % the points/vertices in the deformed cylinder
rad = zeros(a,1);
t = 0;  % index of the current vertex
N = 0;  % number of total vertices
points = true(np,1);
for i = 1:n
    % Modify the refence according to the cylinder axis
    if i > 1
        D = cross(Axe(i,:),U);
        if norm(D) > 0
            a = acos(Axe(i,:)*U);
            R = rotation_matrix(D,-a);
            v = R*V;
            w = R*W;
        else
            v = V;
            w = W;
        end
        D = D/norm(D);
    else
        v = V;
        w = W;
    end
    
    % Define the cylindrical segments
    [d,X,h,Y] = distances_to_line(P(points,:),Axe(i,:),Sta(i,:));
    
    if i > 1
        B = cross(Axe(i-1,:),D);
        b = sum(B.*Axe(i,:));
        if b < 0 || b > pi/2
            B = -B;
        end
        [~,~,hb] = distances_to_line(Sta(i,:)+Rad(i-1)*B,Axe(i,:),Sta(i,:));
    end
    
    J = h < Len(i);   % segment contains points up to the cylinder top
    s = S(points);
    points(s(J)) = false;
    if i > 1
        I = h > hb+H/2;   % segment contains points down to the previous cylinder top + half of H
        J = I&J;
        L = Len(i)-hb-H/2;
        h = h(J)-hb-H/2;
    else
        L = Len(i);
        I = h > 0;
        J = I&J;
        h = h(J);
    end
    d = d(J);       % distances of the points to the cylinder axis
    X = X(J,:);     % vectors connecting the points to the axis
    Y = Y(J,:);     % vectors along the axis to the base of "X"-vectors    
    
    % Define the angles for the points
    D = X*[v w];
    ang = atan2(D(:,1),D(:,2))+pi;  % the angles of the points
    m = ceil(L/H)-1;
    m0 = m;
    
    % Define the layers
    for j = 1:m
        J = h >= (j-1)/m*L;
        K = h <= j/m*L;
        J = J&K;
        % Define the sectors for the layer
        if any(J)
            Ang = ang(J);
            dis = d(J);
            y = Y(J,:);
            x = X(J,:);
            Cent = zeros(NA,3);
            I = false(NA,1);
            for k = 1:NA
                J = Ang < k/NA*2*pi;
                K = Ang >= (k-1)/NA*2*pi;
                J = J&K;
                if any(J)
                    if nnz(J) > 1
                        r = mean(dis(J));
                        C = mean(y(J,:));
                        R = mean(x(J,:));
                    else
                        C = y(J,:);
                        r = dis(J);
                        R = x(J,:);
                    end
                    t = t+1;
                    Vert(t,:) = Sta(i,:)+C+R;  % New vertex
                    rad(t) = r;
                    Cent(k,:) = C;
                    I(k) = true;
                else
                    t = t+1;
                end
            end
            if any(I)
                C = Sta(i,:)+mean(Cent(I,:));
            else
                C = Sta(i,:)+j/m*L*Axe(i,:);
            end
            t = t+1;
            Vert(t,:) = C;  % the center vertex of the layer
            rad(t) = 1;
        else
            m0 = m0-1;
        end
    end
    N = N+(NA+1)*m0;
    
end
Vert = Vert(1:N,:);
rad = rad(1:N);


%% Interpolate the missing surface vertices
M = N/(NA+1);  % total number of vertex layers
NoM = nnz(rad == 0);
NoMpre = N;
while NoM > 0 && NoMpre > NoM
    i = 1;
    while i <= N
        if rad(i) == 0
            m = ceil(i/(NA+1)); % the index of the layer under modification
            
            % Determine the "boundary" vertices "a" and "b" for the missing vertices
            a = i-1;
            if a-(m-1)*(NA+1) == 0 || rad(a) == 0
                a = m*(NA+1)-1;
                while rad(a) == 0
                    a = a-1;
                end
            end
            b = i+1;
            while b <= m*(NA+1)-1 && rad(b) == 0
                b = b+1;
            end
            if b == m*(NA+1)
                b = (m-1)*(NA+1)+1;
                while b <= m*(NA+1)-1 && rad(b) == 0
                    b = b+1;
                end
            end
            if b > a
                K = a+1:1:b-1;
                c = b-a-1;
                r = linspace(rad(a),rad(b),c+2);
                r = r(2:end-1)';
            else
                if a+1 == m*(NA+1)
                    c = NA-(a-b)-1;
                    K = (m-1)*(NA+1)+1:1:b-1;
                elseif b < i && b == (m-1)*(NA+1)+1
                    c = NA-(a-b)-1;
                    K = a+1:1:m*(NA+1)-1;
                elseif b < 1
                    K = [a+1:1:m*(NA+1)-1 (m-1)*(NA+1)+1:1:b-1];
                    c = NA-(a-b)-1;
                else
                    K = [a+1:1:m*(NA+1)-1 (m-1)*(NA+1)+1:1:b-1];
                    c = NA-(a-b)-1;
                end
                r = linspace(rad(a),rad(b),c+2);
                r = r(2:end-1)';
            end
            
            % Add the missing vertices
            if m > 1 && m < M  % the "between" or "inner" layers
                U = K+NA+1;
                L = K-NA-1;
                J1 = rad(U) > 0;
                J2 = rad(L) > 0;
                J = J1&J2;
                if any(J)  % if below and above vertices exists, use them to interpolate between point
                    K = K(J);
                    U = U(J);
                    L = L(J);
                    V = Vert(U,:)-Vert(L,:);
                    V = V/2;
                    Vert(K,:) = Vert(L,:)+V;
                    V = mat_vec_subtraction(Vert(K,:),Vert(m*(NA+1),:));
                    L = sqrt(sum(V.*V,2));
                    rad(K) = L;
                    
                elseif c >= NA/4 % more than a quartan missing, seek the below and above vertices for the first missing vertice
                    a = i;
                    A = a+NA+1;
                    t = 1;
                    while A <= N && rad(A) == 0
                        A = A+NA+1;
                        t = t+1;
                    end
                    B = a-NA-1;
                    if A < N && rad(B) > 0
                        V = Vert(A,:)-Vert(B,:);
                        for j = 1:t
                            Vert(a+(j-1)*(NA+1),:) = Vert(B,:)+j/(t+1)*V;
                            W = Vert(a+(j-1)*(NA+1),:)-Vert((m+j-1)*(NA+1),:);
                            rad(a+(j-1)*(NA+1)) = norm(W);
                        end
                    end
                    
                else % Use rotation, "circle"-approximation
                    U = Vert(a,:)-Vert(m*(NA+1),:);
                    U = U/norm(U);
                    V = Vert(b,:)-Vert(m*(NA+1),:);
                    V = V/norm(V);
                    alp = acos(U*V');
                    D = cross(U,V);
                    vert = zeros(c,3);
                    for j = 1:c
                        R = rotation_matrix(D,alp/(c+1)*j);
                        vert(j,:) = r(j)*(R*U')';
                    end
                    Vert(K,:) = mat_vec_subtraction(vert,-Vert(m*(NA+1),:));
                    rad(K) = r;
                end
            elseif m == 1  % bottom layer
                U = K+NA+1;
                J = rad(U) > 0;
                if any(J)
                    K = K(J);
                    U = U(J);
                    A = rad(1:NA) > 0;
                    B = rad(NA+2:2*NA+1) > 0;
                    J = A&B;
                    A = 1:1:NA;
                    A = A(J);
                    B = NA+2:1:2*NA+1;
                    B = B(J);
                    V = Vert(B,:)-Vert(A,:);
                    if nnz(J) > 1
                        V = mean(V);
                    elseif isempty(V)
                        V = [0 0 H];
                    end
                    Vert(K,:) = mat_vec_subtraction(Vert(U,:),V);
                    V = mat_vec_subtraction(Vert(K,:),Vert(NA+1,:));
                    L = sqrt(sum(V.*V,2));
                    rad(K) = L;
                else
                    % Use rotation, "circle"-approximation
                    U = Vert(a,:)-Vert(NA+1,:);
                    U = U/norm(U);
                    V = Vert(b,:)-Vert(NA+1,:);
                    V = V/norm(V);
                    alp = acos(U*V');
                    D = cross(U,V);
                    vert = zeros(c,3);
                    for j = 1:c
                        R = rotation_matrix(D,alp/(c+1)*j);
                        vert(j,:) = r(j)*(R*U')';
                    end
                    Vert(K,:) = mat_vec_subtraction(vert,-Vert(NA+1,:));
                    rad(K) = r;
                end
            else  % top layer
                L = K-NA-1;
                J = rad(K) > 0;
                L = L(J);
                K = K(J);
                J = rad((m-1)*(NA+1)+1:m*(NA+1)-1) > 0;
                A = (m-1)*(NA+1)+1:1:m*(NA+1)-1;
                A = A(J);
                B = A-NA-1;
                V = Vert(B,:)-Vert(A,:);
                if nnz(J) > 1
                    V = mean(V);
                end
                Vert(K,:) = mat_vec_subtraction(Vert(L,:),V);
                V = mat_vec_subtraction(Vert(K,:),Vert(m*(NA+1),:));
                L = sqrt(sum(V.*V,2));
                rad(K) = L;
            end
            i = i+1;
        else
            i = i+1;
        end
    end
    
    % Update the information of number of missing vertices
    NoMpre = NoM;
    NoM = nnz(rad == 0);
end


%% Make the top layer complete
K = rad(N-NA:N) == 0;
while any(K)
    i = N-NA;
    while i <= N
        if rad(i) == 0
            m = ceil(i/(NA+1)); % the index of the layer under modification
            
            % Determine the boundary vertices "a" and "b" for the missing vertices
            a = i-1;
            if a-(m-1)*(NA+1) == 0 || rad(a) == 0
                a = m*(NA+1)-1;
                while rad(a) == 0
                    a = a-1;
                end
            end
            b = i+1;
            while b <= N && rad(b) == 0
                b = b+1;
            end
            if b == m*(NA+1)
                b = (m-1)*(NA+1)+1;
            end
            if b > N
                b = N;
            end
            
            if b > a
                K = a+1:1:b-1;
                c = b-a-1;
                r = linspace(rad(a),rad(b),c+2);
                r = r(2:end-1)';
            else
                if a+1 == m*(NA+1)
                    c = NA-(a-b)-1;
                    K = (m-1)*(NA+1)+1:1:b-1;
                elseif b < i && b == (m-1)*(NA+1)+1
                    c = NA-(a-b)-1;
                    K = a+1:1:m*(NA+1)-1;
                elseif b < 1
                    K = [a+1:1:m*(NA+1)-1 (m-1)*(NA+1)+1:1:b-1];
                    c = NA-(a-b)-1;
                else
                    K = [a+1:1:m*(NA+1)-1 (m-1)*(NA+1)+1:1:b-1];
                    c = NA-(a-b)-1;
                end
                r = linspace(rad(a),rad(b),c+2);
                r = r(2:end-1)';
            end
                        
            U = K-NA-1;
            J = rad(U) > 0;
            if any(J)
                K = K(J);
                U = U(J);
                JU = N-NA:1:N-1;
                JL = JU-NA-1;
                A = rad(JL) > 0;
                B = rad(JU) > 0;
                J = A&B;
                if any(J)
                    A = JL(J);
                    B = JU(J);
                    V = Vert(B,:)-Vert(A,:);
                    if nnz(J) > 1
                        V = mean(V);
                    end
                    Vert(K,:) = mat_vec_subtraction(Vert(U,:),-V);
                    V = mat_vec_subtraction(Vert(K,:),Vert(N,:));
                    L = sqrt(sum(V.*V,2));
                    rad(K) = L;
                else
                    % Use rotation, "circle"-approximation
                    U = Vert(a,:)-Vert(N,:);
                    U = U/norm(U);
                    if b ~= a
                        V = Vert(b,:)-Vert(N,:);
                        V = V/norm(V);
                    else
                        b = a-1;
                        V = Vert(b,:)-Vert(N,:);
                        V = V/norm(V);
                    end
                    alp = acos(U*V');
                    D = cross(U,V);
                    vert = zeros(c,3);
                    for j = 1:c
                        R = rotation_matrix(D,alp/(c+1)*j);
                        vert(j,:) = r(j)*(R*U')';
                    end
                    Vert(K,:) = mat_vec_subtraction(vert,-Vert(N,:));
                    rad(K) = r;
                end
            else
                % Use rotation, "circle"-approximation
                U = Vert(a,:)-Vert(N,:);
                U = U/norm(U);
                V = Vert(b,:)-Vert(N,:);
                V = V/norm(V);
                alp = acos(U*V');
                D = cross(U,V);
                vert = zeros(c,3);
                for j = 1:c
                    R = rotation_matrix(D,alp/(c+1)*j);
                    vert(j,:) = r(j)*(R*U')';
                end
                Vert(K,:) = mat_vec_subtraction(vert,-Vert(N,:));
                rad(K) = r;
            end
            if rad(i) > 0
                i = i+1;
            end
        else
            i = i+1;
        end
    end
    K = rad(N-NA:N) == 0;
end


%% Interpolate if still missing vertices
NoM = nnz(rad == 0);
NoMpre = N;
while NoM > 0 && NoMpre > NoM
    i = 1;
    while i <= N
        if rad(i) == 0
            m = ceil(i/(NA+1)); % the index of the layer under modification
            
            % Determine the boundary vertices "a" and "b" for the missing vertices
            a = i-1;
            if a-(m-1)*(NA+1) == 0 || rad(a) == 0
                a = m*(NA+1)-1;
                while rad(a) == 0
                    a = a-1;
                end
            end
            b = i+1;
            while b <= m*(NA+1)-1 && rad(b) == 0
                b = b+1;
            end
            if b == m*(NA+1)
                b = (m-1)*(NA+1)+1;
                while b <= m*(NA+1)-1 && rad(b) == 0
                    b = b+1;
                end
            end
            
            if b > a
                K = a+1:1:b-1;
                c = b-a-1;
                r = linspace(rad(a),rad(b),c+2);
                r = r(2:end-1)';
            else
                if a+1 == m*(NA+1)
                    c = NA-(a-b)-1;
                    K = (m-1)*(NA+1)+1:1:b-1;
                elseif b < i && b == (m-1)*(NA+1)+1
                    c = NA-(a-b)-1;
                    K = a+1:1:m*(NA+1)-1;
                elseif b < 1
                    K = [a+1:1:m*(NA+1)-1 (m-1)*(NA+1)+1:1:b-1];
                    c = NA-(a-b)-1;
                else
                    K = [a+1:1:m*(NA+1)-1 (m-1)*(NA+1)+1:1:b-1];
                    c = NA-(a-b)-1;
                end
                r = linspace(rad(a),rad(b),c+2);
                r = r(2:end-1)';
            end
            
            % Add the missing vertices
            if m > 1 && m < M  % the "between" or "inner" layers
                U = K+NA+1;
                L = K-NA-1;
                J1 = rad(U) > 0;
                J2 = rad(L) > 0;
                J = J1&J2;
                if any(J)
                    K = K(J);
                    U = U(J);
                    L = L(J);
                    V = Vert(U,:)-Vert(L,:);
                    V = V/2;
                    Vert(K,:) = Vert(L,:)+V;
                    V = mat_vec_subtraction(Vert(K,:),Vert(m*(NA+1),:));
                    L = sqrt(sum(V.*V,2));
                    rad(K) = L;
                elseif c >= NA/4
                    a = i;
                    A = a+NA+1;
                    t = 1;
                    while A <= N && rad(A) == 0
                        A = A+NA+1;
                        t = t+1;
                    end
                    B = a-NA-1;
                    if A < N && rad(B) > 0
                        V = Vert(A,:)-Vert(B,:);
                        for j = 1:t
                            Vert(a+(j-1)*(NA+1),:) = Vert(B,:)+j/(t+1)*V;
                            W = Vert(a+(j-1)*(NA+1),:)-Vert((m+j-1)*(NA+1),:);
                            rad(a+(j-1)*(NA+1)) = norm(W);
                        end
                    end
                    
                else
                    % Use rotation, "circle"-approximation
                    U = Vert(a,:)-Vert(m*(NA+1),:);
                    U = U/norm(U);
                    V = Vert(b,:)-Vert(m*(NA+1),:);
                    V = V/norm(V);
                    alp = acos(U*V');
                    D = cross(U,V);
                    vert = zeros(c,3);
                    for j = 1:c
                        R = rotation_matrix(D,alp/(c+1)*j);
                        vert(j,:) = r(j)*(R*U')';
                    end
                    Vert(K,:) = mat_vec_subtraction(vert,-Vert(m*(NA+1),:));
                    rad(K) = r;
                end
            elseif m == 1  % bottom layer
                U = K+NA+1;
                J = rad(U) > 0;
                if any(J)
                    K = K(J);
                    U = U(J);
                    A = rad(1:NA) > 0;
                    B = rad(NA+2:2*NA+1) > 0;
                    J = A&B;
                    A = 1:1:NA;
                    A = A(J);
                    B = NA+2:1:2*NA+1;
                    B = B(J);
                    V = Vert(B,:)-Vert(A,:);
                    if nnz(J) > 1
                        V = mean(V);
                    end
                    Vert(K,:) = mat_vec_subtraction(Vert(U,:),V);
                    V = mat_vec_subtraction(Vert(K,:),Vert(NA+1,:));
                    L = sqrt(sum(V.*V,2));
                    rad(K) = L;
                else
                    % Use rotation, "circle"-approximation
                    U = Vert(a,:)-Vert(NA+1,:);
                    U = U/norm(U);
                    V = Vert(b,:)-Vert(NA+1,:);
                    V = V/norm(V);
                    alp = acos(U*V');
                    D = cross(U,V);
                    vert = zeros(c,3);
                    for j = 1:c
                        R = rotation_matrix(D,alp/(c+1)*j);
                        vert(j,:) = r(j)*(R*U')';
                    end
                    Vert(K,:) = mat_vec_subtraction(vert,-Vert(NA+1,:));
                    rad(K) = r;
                end
            else  % top layer
                L = K-NA-1;
                J = rad(K) > 0;
                L = L(J);
                K = K(J);
                J = rad((m-1)*(NA+1)+1:m*(NA+1)-1) > 0;
                A = (m-1)*(NA+1)+1:1:m*(NA+1)-1;
                A = A(J);
                B = A-NA-1;
                V = Vert(B,:)-Vert(A,:);
                if nnz(J) > 1
                    V = mean(V);
                end
                Vert(K,:) = mat_vec_subtraction(Vert(L,:),V);
                V = mat_vec_subtraction(Vert(K,:),Vert(m*(NA+1),:));
                L = sqrt(sum(V.*V,2));
                rad(K) = L;
            end
            i = i+1;
        else
            i = i+1;
        end
    end
    
    % Update the information of number of missing vertices
    NoMpre = NoM;
    NoM = nnz(rad == 0);
end

%% Problems probably with the top part, make the reconstruction shorter
if NoM > 0
    i = 1;
    while rad(i) > 0
        i = i+1;
    end
    M = ceil(i/(NA+1))-1;
end


%% Try to compensate bottom and top vertices
V = mat_vec_subtraction(Vert(1:NA+1,:),H/2*Axe(1,:));
Vert(1:NA+1,:) = V;
V = mat_vec_subtraction(Vert(end-NA:end,:),-H/2*Axe(n,:));
Vert(end-NA:end,:) = V;


%% Diameter data
m = ceil(N/(NA+1));
DiamSummary = zeros(m,4);
h = 0;
for i = 1:m
    V = Vert((i-1)*(NA+1)+NA/2+1:i*(NA+1)-1,:)-Vert((i-1)*(NA+1)+1:(i-1)*(NA+1)+NA/2,:);
    d = sqrt(sum(V.*V,2));
    if i > 1
        h = h+norm(Vert(i*(NA+1),:)-Vert((i-1)*(NA+1),:));
    end
    DiamSummary(i,:) = [min(d) mean(d) max(d) h];
end
DiamSummary = 100*DiamSummary;


%% Check if "blow-up" and shorten
i = 3;
D = DiamSummary(:,2);
while i <= m && D(i) < 1.25*D(i-1) && D(i) < 1.25*D(i-2)
    i = i+1;
end
if i <= m
    M = i-2;
    Vert = Vert(1:(NA+1)*M,:);
    DiamSummary = DiamSummary(1:M,:);
    disp('Shortening the triangulated surface')
end


%% Define the triangles which build up the surface and its volume
m = M-1;
Tria = zeros(NA*m,4);
if m > 1
    Bot = zeros(NA,3);
    Top = zeros(NA,3);
else
    Bot = zeros(0,3);
    Top = zeros(0,3);
end
for i = 1:m
    if i == 1
        for j = 1:NA
            if j < NA
                Bot(j,:) = [NA+1 j j+1];
            else
                Bot(j,:) = [NA+1 j 1];
            end
        end
    elseif i == m
        a = (m+1)*(NA+1);
        b = m*(NA+1);
        for j = 1:NA
            if j < NA
                Top(j,:) = [a b+j b+j+1];
            else
                Top(j,:) = [a b+j b+1];
            end
        end
    end
    a = i*(NA+1);
    b = (i-1)*(NA+1);
    for j = 1:NA
        if j < NA
            Tria((i-1)*NA+j,:) = [b+j b+j+1 a+j a+j+1];
        else
            Tria((i-1)*NA+j,:) = [b+j b+1 a+j a+1];
        end
    end
end


%% Determine the volume
% The volume is calculated by the formula derived from the divergence
% theorem, where the volume is certain sum over the facets.

if m > 1
    % Side
    U = Vert(Tria(:,3),:)-Vert(Tria(:,1),:);
    V = Vert(Tria(:,4),:)-Vert(Tria(:,1),:);
    W = Vert(Tria(:,2),:)-Vert(Tria(:,1),:);
    C1 = Vert(Tria(:,1),:)+0.25*V+0.25*U;
    C2 = Vert(Tria(:,1),:)+0.25*V+0.25*W;
    N1 = cross(U,V);
    N2 = cross(V,W);
    A1 = 0.5*sqrt(sum(N1.*N1,2));
    A2 = 0.5*sqrt(sum(N2.*N2,2));
    N1 = 0.5*[N1(:,1)./A1 N1(:,2)./A1 N1(:,3)./A1];
    N2 = 0.5*[N2(:,1)./A2 N2(:,2)./A2 N2(:,3)./A2];
    
    % Bottom
    V = Vert(Bot(:,3),:)-Vert(Bot(:,1),:);
    U = Vert(Bot(:,2),:)-Vert(Bot(:,1),:);
    C3 = Vert(Bot(:,1),:)+0.25*V+0.25*U;
    N3 = cross(U,V);
    A3 = 0.5*sqrt(sum(N3.*N3,2));
    N3 = 0.5*[N3(:,1)./A3 N3(:,2)./A3 N3(:,3)./A3];
    
    % Top
    V = Vert(Top(:,3),:)-Vert(Top(:,1),:);
    U = Vert(Top(:,2),:)-Vert(Top(:,1),:);
    C4 = Vert(Top(:,1),:)+0.25*V+0.25*U;
    N4 = cross(V,U);
    A4 = 0.5*sqrt(sum(N4.*N4,2));
    N4 = 0.5*[N4(:,1)./A4 N4(:,2)./A4 N4(:,3)./A4];
    
    % Compute the volume (in liters) and round up
    Vtot = sum(A1.*sum(C1.*N1,2))+sum(A2.*sum(C2.*N2,2))+...
        sum(A3.*sum(C3.*N3,2))+sum(A4.*sum(C4.*N4,2));
    Vtot = 1000*Vtot/3;
    if Vtot > 100
        Vtot = round(Vtot);
    elseif Vtot > 10
        Vtot = round(10*Vtot)/10;
    else
        Vtot = round(100*Vtot)/100;
    end
    
    Facets = [Tria(:,1) Tria(:,3) Tria(:,4); Tria(:,1) Tria(:,4) Tria(:,2); ...
        Bot(:,1) Bot(:,2) Bot(:,3); Top(:,1) Top(:,2) Top(:,3)];
    fvd = ones(size(Facets,1),1);
    
%     figure(4)
%     plot3(P(:,1),P(:,2),P(:,3),'.r','Markersize',7)
%     patch('Vertices',Vert,'Faces',Facets,'FaceVertexCData',fvd,'FaceColor','flat')
%     axis equal
%     alpha(0.95)
%     pause(0.3)
    
else
    Vtot = 0;
end