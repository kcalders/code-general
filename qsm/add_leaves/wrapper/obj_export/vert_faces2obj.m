function vert_faces2obj(Vertices,Faces,d,file)

    % Set precision formatter.
    if length(d) > 1
        ft = ['%' num2str(d(2)) '.' num2str(d(1)) 'f'];
    else
        ft = ['%.' num2str(d(1)) 'f'];
    end

    % Flag to close file and the end.
    closefile = false;

    % Check if filename and not file stream.
    if ischar(file)

        % Open file stream with filename.
        fid = fopen(file,'w');
        % Set file to close at the end.
        closefile = true;

    else
        % Otherwise a file stream is given as input.
        fid = file;
    end

    NVertex = size(Vertices,1);

    % Print vertices.
    for iVertex = 1:NVertex
        
        % Index <leaves> and move to origin.
        fprintf(fid,['v ' ft ' ' ft ' ' ft '\n'],Vertices(iVertex,:));
        
    end

    NFace = size(Faces,1);

    fCellFaces = false;

    if iscell(Faces)
        fCellFaces = true;
    end

    % Print faces.
    for iFace = 1:NFace

        if fCellFaces

            % Vertices in face
            NVertex = length(Faces{iFace});

            ft = repmat('%d ',1,NVertex);
            ft = ft(1:end);

            fprintf(fid,['f ' ft '\n'],Faces{iFace});

        else
            if iFace == 1
                NVertex = size(Faces,2);

                ft = repmat('%d ',1,NVertex);
                ft = ft(1:end);
            end

            fprintf(fid,['f ' ft '\n'],Faces(iFace,:));
        end

    end

end