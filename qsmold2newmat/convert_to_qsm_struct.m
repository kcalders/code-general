function QSMs = convert_to_qsm_struct(folder,ID,Species)

% Given the folder name as a string "folder", the code reconstruct
% QSMs structure array similar to TreeQSM (version >= 2.3.0) from
% the outputs of the older TreeQSM versions.
% Optional inputs "ID" and "Species", which are strings of tree id:s and
% their species labels, add additional field "species" to the output 
% structure array QSMs.

tic
D = dir(folder);
n = max(size(D))
QSMs = struct('cylinder',{},'branch',{},'treedata',{},'name',{},'species',{});

j = 0;
inputs.disp = 0;
inputs.plot = 0;
c0 = 1;

for i = 1:n
    if D(i).bytes > 1e5
        % Load the tree files:
        load([D(i).folder,'/',D(i).name])
        
        % Reconstruct the cylinder structure:
        cylinder.radius = single(Rad);
        cylinder.length = single(Len);
        cylinder.axis = single(Axe);
        cylinder.start = single(Sta);
        m = length(Rad);
        if m <= 2^16
            cylinder.parent = uint16(CPar);
            cylinder.extension = uint16(CExt);
        else
            cylinder.parent = uint32(CPar);
            cylinder.extension = uint32(CExt);   
        end
        cylinder.branch = uint16(BoC(:,1));
        cylinder.BranchOrder = uint8(BoC(:,2));
        if max(BoC(:,3)) <= 2^8
            cylinder.PositionInBranch = uint8(BoC(:,3));
        else
            cylinder.PositionInBranch = uint16(BoC(:,3));
        end
        cylinder.added = Added;
        
        % Reconstruct branch and treedata structures:
        branch = branches2(cylinder);
        treedata = tree_data2(cylinder,branch,inputs);
        
        % Add the structures to the output:
        j = j+1;
        QSMs(j).cylinder = cylinder;
        QSMs(j).branch = branch;
        QSMs(j).treedata = treedata;
        QSMs(j).name = D(i).name;
        
        if nargin > 1
            % Add the species label:
            % Search from the file name the part corresponding to the tree
            % IDs in the input string array "ID". File names are now
            % assumed to be in the form "xxxxx_xxxx_treeid-xxx-xxx..."
            name = D(i).name;
            b = 1;
            while ~strcmp(name(b),'-')
                b = b+1;
            end
            b = b-1;
            %QSMs(j).name = name(1:k);
            a = b;
            while ~strcmp(name(a),'_')
                a = a-1;
            end
            a = a+1;
            
            % Search the same treeid from the input "ID" to find the
            % correct species label:
            c = c0;
            treeid = name(a:b);
            while ~strcmp(ID(c,:),treeid)
                c = c+1;
            end
            if c == c0
                c0 = c0+1;
            end
            % Add the species label:
            QSMs(j).species = Species(c);
        end
        
        if floor(j/10) == ceil(j/10)
            disp(j)
            toc
        end
    end
end