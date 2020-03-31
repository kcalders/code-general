function [Segs,SPar,SChi,SDir] = correct_segments(P,Cen,Segs,SPar,SChi,SDir)

% ---------------------------------------------------------------------
% CORRECT_SEGMENTS.M        Combines and removes segments.
%
% Version 1.0
% Author        Pasi Raumonen
% Created       14 June 2013
%
% This software cannot be distributed without the permission of P. Raumonen
% and it is only for non-commercial use. These restrictions holds also for
% the modules or subprograms of the software:
% DISTANCES_TO_LINE.M
% ---------------------------------------------------------------------

% Combine segments to form new segments. If a small segment is along a
% larger segment as a child segment, then remove the smaller one. If a
% child segment continues its parent segment, then combine theM to a single
% longer segment.


%% Remove small child segments
ns = size(Segs,1);
I = true(ns,1);
for i = 1:ns
    C = SChi{i};
    if ~isempty(C)
        n = length(C);
        S = Segs{i};
        s = size(S,1);
        for j = 1:n
            nl = SPar(C(j),2);
            if nl > 10
                a = nl-10;
            else
                a = 1;
            end
            if s-nl > 10
                b = nl+10;
            else
                b = s;
            end
            B = S{b};
            if length(B) > 1
                B = mean(P(Cen(B),:));
            else
                B = P(Cen(B),:);
            end
            A = S{a};
            if length(A) > 1
                A = mean(P(Cen(A),:));
            else
                A = P(Cen(A),:);
            end
            V = B-A;
            V = V/norm(V);
            Q = vertcat(S{a:b});
            R = Segs{C(j)};
            R = vertcat(R{:});
            R = [Q; R];
            dR = max(distances_to_line(P(Cen(R),:),V,A));
            dQ = max(distances_to_line(P(Cen(Q),:),V,A));
            if (dR-dQ < 0.02) || (dR/dQ < 1.2 && dR-dQ < 0.06)
                K = SChi{C(j)};
                c = length(K);
                if isempty(K)
                    % Remove, no child segments
                    I(C(j)) = false; 
                    Segs{C(j)} = zeros(0,1);
                    SPar(C(j),:) = zeros(1,2);
                    SChi{i} = setdiff(C,C(j));
                else
                    L = SChi(K);
                    L = vertcat(L{:}); % child child segments
                    if isempty(L)
                        J = false(c,1);
                        for k = 1:c
                            r = Segs{K(k)};
                            r = vertcat(r{:});
                            r = [r; R];
                            dr = max(distances_to_line(P(Cen(r),:),V,A));
                            if (dr-dQ < 0.02) || (dr/dQ < 1.2 && dr-dQ < 0.06)
                                J(k) = true;
                            end
                        end
                        if all(J)
                            % Remove
                            K = [K; C(j)];
                            c = length(K);
                            Segs(K) = cell(c,1);
                            I(K) = false; 
                            SPar(K,:) = zeros(c,2);
                            SChi{i} = setdiff(C,C(j));
                            SChi(K) = cell(c,1);
                        end
                    end
                end
            end
        end
    end
end
Segs = Segs(I);
Ind = (1:1:ns)';
n = nnz(I);
J = (1:1:n)';
Ind(I) = J;
Ind(~I) = 0;
SPar = SPar(I,:);
SDir = SDir(I,:);
J = SPar(:,1) > 0;
SPar(J,1) = Ind(SPar(J,1));
% Modify SChi
for i = 1:ns
    if I(i)
        C = SChi{i};
        if ~isempty(C)
            C = nonzeros(Ind(C));
            SChi{i} = C;
        end
    end
end
SChi = SChi(I);

ns = size(Segs,1);

%% Combine segments
ns0 = ns+1;
while ns < ns0
    E = true(ns,1);
    for i = 2:ns
        C = SChi{i};
        S = Segs{i};
        n = size(S,1);
        s = SPar(C,2);
        I = s >= n-1;
        C = C(I);
        if ~isempty(C) && ~isempty(S)
            if n > 5
                S0 = S{n-5};
            else
                S0 = S{1};
            end
            if length(S0) > 1
                Q0 = mean(P(Cen(S0),:));
            else
                Q0 = P(Cen(S0),:);
            end
            S = S{end};
            if length(S) > 1
                Q = mean(P(Cen(S),:));
            else
                Q = P(Cen(S),:);
            end
            V = Q-Q0;
            V = V/norm(V);
            nc = length(C);
            J = false(nc,1);
            q = zeros(nc,3);
            for j = 1:nc
                s = Segs{C(j)};
                if ~isempty(s)
                    if length(s) > 5
                        S0 = s{6};
                    else
                        S0 = s{end};
                    end
                    if length(S0) > 1
                        Q = mean(P(Cen(S0),:));
                    else
                        Q = P(Cen(S0),:);
                    end
                    s = s{1};
                    J(j) = true;
                    if length(s) > 1
                        Q0 = mean(P(Cen(s),:));
                    else
                        Q0 = P(Cen(s),:);
                    end
                    W = Q-Q0;
                    W = W/norm(W);
                    q(j,:) = W;
                end
            end
            if any(J)
                I = q*V';
                [m,k] = max(I);
                C = C(k);
                
                S = Segs{i};
                s = Segs{C};
                Segs{i} = [S; s];
                Segs{C} = zeros(0,1);
                c = SChi{i};
                s = SChi{C};
                SPar(s,:) = [i*ones(length(s),1) size(S,1)+SPar(s,2)];
                c = union(c,s);
                c = setdiff(c,C);
                SChi{i} = c;
                E(C) = false;
            end
        end
    end
    
    Segs = Segs(E);
    Ind = (1:1:ns)';
    n = nnz(E);
    I = (1:1:n)';
    Ind(E) = I;
    SPar = SPar(E,:);
    SDir = SDir(E,:);
    J = SPar(:,1) > 0;
    SPar(J,1) = Ind(SPar(J,1));
    for i = 1:ns
        if E(i)
            C = SChi{i};
            if ~isempty(C)
                C = Ind(C);
                SChi{i} = C;
            end
        end
    end
    SChi = SChi(E);
    ns0 = ns;
    ns = nnz(E);
end
SPar = SPar(:,1);

str = ['    ',num2str(ns),' segments after correction'];
disp(str)