function N = LeafModel2obj(Leaves,d,file,OriginOffset,Filter)

    fOriginOverride = false;

    if nargin > 3 && not(isempty(OriginOffset))
        fOriginOverride = true;
    end
    
    if nargin <= 4
        Filter = true(Leaves.leaf_count,1);
    end

    fNgon = true;

    if isempty(Leaves.base_ngons)
        fNgon = false;
    end

    % Number of leaves in model.
    NLeaf = Leaves.leaf_count;
    
    % Number of vertices in base.
    NBaseVert = size(Leaves.base_vertices,1);
    
    % Matrix of vertices, row == vertex.
    Vertices = zeros(nnz(Filter)*NBaseVert,3);
    
    jLeaf = 0;
    
    for iLeaf = 1:NLeaf
        
        if Filter(iLeaf)
            jLeaf = jLeaf + 1;
        else
            continue;
        end

        % Get single leaf parameters.
        origin = Leaves.leaf_start_point(iLeaf,:);
        scale  = Leaves.leaf_scale(iLeaf,:);
        dir    = Leaves.leaf_direction(iLeaf,:);
        normal = Leaves.leaf_normal(iLeaf,:);
        
        % Compute vertices by transforming leaf base.
        vert  = Leaves.base_vertices;

        % Scaling.
        vert = bsxfun(@times,vert,scale);

        % Coordinate system.
        E = [cross(normal,dir); dir; normal];

        % Rotation.
        vert = vert*E;

        % Translation.
        vert = bsxfun(@plus,vert,origin);

        % Store resulting vertices.
        Vertices((jLeaf-1)*NBaseVert+1:jLeaf*NBaseVert,:) = vert;
        
    end
    
    NLeaf = jLeaf;

    % Using ngons
    if fNgon

        fEqualNgons = true;

        try
            BaseNgons = vertcat(Leaves.base_ngons{:});
        catch
            fEqualNgons = false;
        end

        % All ngons have the same number of vertices.
        if fEqualNgons

            % Number of triangles in base.
            NNgon = size(BaseNgons,1);

            % Indices of base triangle face vertices.
            ngons = repmat(BaseNgons,NLeaf,1);

            add = repmat(0:1:NLeaf-1,NNgon,1);
            add = add(:)*NBaseVert;

            Faces = bsxfun(@plus,ngons,add);

        % Ngon vertex count varies.
        else

            % Base ngons as cell array.
            BaseNgons = Leaves.base_ngons;

            % Number of ngons in base.
            NNgon = numel(BaseNgons);

            % Cell to store faces.
            Faces = cell(NLeaf*NNgon,1);

            for iLeaf = 1:NLeaf

                for iNgon = 1:NNgon

                    % Offset index by previus leaf count,
                    % and store to cell array.
                    Faces{(iLeaf-1)*NNgon+iNgon} = BaseNgons{iNgon} ...
                                                 + (iLeaf-1)*NNgon;
                    %-

                end

            end
            
        end

    else
        
        % Number of triangles in base.
        NTri = size(Leaves.base_triangles,1);

        % Indices of base triangle face vertices.
        tris = repmat(Leaves.base_triangles,NLeaf,1);

        % Offset vector to account for previous leaves.
        add = repmat(0:1:NLeaf-1,NTri,1);
        add = add(:)*NBaseVert;

        % Store offset face indices to matrix.
        Faces = bsxfun(@plus,tris,add);

    end

    % Override origin if given.
    if fOriginOverride
        Vertices = bsxfun(@minus,Vertices,OriginOffset);
    end

    % Write resulting vertices and faces to file.
    vert_faces2obj(Vertices,Faces,d,file);
    
    % Return number of written vertices.
    N = size(Vertices,1);

end