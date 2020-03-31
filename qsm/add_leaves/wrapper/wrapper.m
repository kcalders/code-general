projectfolder = 'qsm-fanni/src/';

addpath(projectfolder);
addpath([projectfolder 'classes']);
addpath('distribution_functions');

outputfolder = '';

VerStr = '1.0';

%% Load models to cell array.

inputfolder = 'path/to/modeldatafiles/';

files = dir([inputfile 'Modeldata*']);

NTree = size(files);

QSMs = cell(NTree,1);

for iTree = 1:NTree

    load([inputfolder files(iTree).name]);

    ModelData = {CylData,BranchData,TreeData};

    QSMs{iTree} = ModelData;
    
    clear CylData BranchData TreeData;

end

%%

% Vertices and faces of the leaf basis geometry.
vertices = [0 0 0; -0.30834 0.66549 0; 0 1 0; 0.30834 0.66549 0];
tris = [1 2 3; 1 3 4];
extraparams = {{[1 2 3 4]}};

% Leaf models.
LeafModels =  cell(NTree,1);

% Resulting leaf counts.
LeafCounts = zeros(NTree,4);

% Computation times.
Times      = zeros(NTree,1);

for iTree = 1:NTree

    % Create an object that holds QSM cylinder properties.
%    QSM = QSMBCylindrical(QSMs{iTree});
	QSM = QSMBCylindrical(QSMs{1});

    % Initialize leaf model:
    % 1) Vertices,
    % 2) Triangles of the faces,
    % 3) Ngons of the faces, if any, and
    % 4) Number of elements to initialize for performance.
    Leaves = LeafModelTriangle(vertices, tris, extraparams{:}, 1e5);
    
    % Leaf area to insert:
    % 1) Area amount after which is reached, it is ok to stop.
    % 2) Area amount to initially generate.
    LeafArea = [250 300];

    % Parameters to fine tune leaf insertion.
    LeafParam = {
                    ...% Leaf area density distribution function.
                    'AreaFunction',@ladd_6,...
                    ...% Limits for leaf length.
                    'SizeFunctionParameters', {[0.25 0.30]},...
                    ...% Disable printing.
                    'Verbose',false
                };
    %-
    
    tic;

    [Leaves, NAccepted] = qsm_fanni(QSM,... % Tree model
                                    Leaves,... % Leaf basis model
                                    LeafArea,... % Target area
                                    'Seed',iTree,... % Set seed for reproductability (opt.)
                                    LeafParam{:}); % Fine tuning.
    %-

    Time = toc;
    
    % Remove empty rows from matrices that might have been initialized.
    Leaves.trim_slack();
    
    % Store results.
    LeafModels{iTree} = Leaves;
    LeafCounts(iTree) = NAccepted;
    Times(iTree) = Time;
end

%% Write leaf models to file as OBJ-files.

addpath('obj_export');

[status, ~,~] = mkdir([outputfolder 'export/'],VerStr);
folder = [outputfolder 'export/' VerStr '/'];

for iTree = 1:NTree
    
    leaffile = [folder 'LeafModel_' num2str(iTree) ...
                '_NLeaf_' num2str(LeafCounts(iTree)) ...
                '.obj'];
    %-

    %
    LeafModel2obj(LeafModels{iTree},4,leaffile);

end

exit