function NP = ladd_6(NBlock,BlocksParameters,varargin)
% Two parameter LADD.

    if nargin > 2
        y0 = varargin{1};
    else
        y0 = 0.2;
    end
    
    if nargin > 3
        y4 = varargin{2};
    else
        y4 = 0.7;
    end

    % Compute leaf distribution parameters for each cylinder.
    NP = zeros(NBlock,1);

    bo = BlocksParameters.branch_order;
    relpos = BlocksParameters.relative_position;
    relheight = BlocksParameters.relative_height;

    HeightPoly = polyfit([0,1],[y0,1],1);

    CutOffPoly = polyfit([0 4],[0.95 y4],1);

    for iBlock = 1:NBlock

        CutOff = polyval(CutOffPoly,min(bo(iBlock),4));

        if relpos(iBlock) > CutOff
            
            PolCoeff = polyfit([CutOff 1],[0 1],1);

            P = polyval(PolCoeff,relpos(iBlock));

            NP(iBlock) = P*polyval(HeightPoly,relheight(iBlock));

        end

    end
    
    NP = NP/sum(NP);

end