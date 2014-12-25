function Plot(obj)
cartesian = [reshape(obj.Molecule_AtomicNumbers(), [], 1), ...
    obj.Molecule_Geometry() .* 0.529177249];

distMat = dist(cartesian(:, 2:end)');
numAtoms = size(cartesian, 1);

maxBonds = zeros(numAtoms, 1);
for iAtom = 1:numAtoms
    maxBonds(iAtom) = MaxNumBonds(cartesian(iAtom, 1));
end

numBondsMat = zeros(numAtoms);
for iter = [1 100]
    for iAtom = 1:numAtoms
        [~, tempOrder] = sort(distMat(iAtom, :));
        for jAtom = tempOrder(2:end)
            if(sum(numBondsMat(iAtom, :)) < maxBonds(iAtom) ...
                    && sum(numBondsMat(jAtom,:)) < maxBonds(jAtom))
                maxBetween = min(maxBonds(iAtom), maxBonds(jAtom));
                numBonds = min(min(NumBonds(cartesian, iAtom, jAtom), maxBetween), iter);
                numBondsMat(iAtom, jAtom) = numBonds;
                numBondsMat(jAtom, iAtom) = numBonds;
            end
        end
    end
end

bondLines = {};
for iAtom = 1:numAtoms
    iPos = cartesian(iAtom, 2:end);
    for jAtom = iAtom+1:numAtoms
        jPos = cartesian(jAtom, 2:end);
        numBonds = numBondsMat(iAtom, jAtom);
        vector = iPos - jPos;
        normVec = ones(1, 3);
        normVec(3) = -(vector(1) + vector(2))/vector(3);
        normVec = normVec ./ norm(normVec);
        ijVecs = [iPos; jPos];
        normVecs = 0.05 .* [normVec; normVec];
        switch(numBonds)
            case(1) % single bond
                bondLines{end+1} = ijVecs;
            case(2) % double bond
                bondLines{end+1} = ijVecs - normVecs;
                bondLines{end+1} = ijVecs + normVecs;
            case(3)% triple bond
                bondLines{end+1} = ijVecs;
                bondLines{end+1} = ijVecs - normVecs;
                bondLines{end+1} = ijVecs + normVecs;
        end
    end
end

figure();
set(gcf,'color','w');
axis vis3d;
axis off;
hold;
for i = 1:length(bondLines)
    line(bondLines{i}(:, 1), bondLines{i}(:, 2), bondLines{i}(:, 3), 'LineWidth', 4);
end

allShapes = 4000 .* [ ...
    1             1 ...
    3 3 2 2 2 2 2 2 ...
    4 4 3 3 3 3 3 3 ...
    5 5];
allColors = 0.12 .* [ ...
    [7 7 7];                                                       [8 7 7]; ...
    [0 8 8]; [8 0 8]; [8 8 0]; [5 5 5]; [0 0 8]; [8 0 0]; [0 8 0]; [6 5 5];
    [0 6 6]; [6 0 6]; [6 6 0]; [3 3 3]; [0 0 6]; [6 0 0]; [0 6 0]; [4 3 3];
    [0 4 4]; [4 0 4];];

scatter3(cartesian(:, 2), cartesian(:, 3), cartesian(:, 4), ...
    allShapes(cartesian(:, 1)), allColors(cartesian(:, 1), :), 'fill');

end


function maxNumBonds = MaxNumBonds(atomicNum)
[~, col, rowMax] = FromAtomicNumber(atomicNum);
maxNumBonds = min(col, rowMax-col);
end

function numBonds = NumBonds(cartesian, indAtom1, indAtom2)

distance = norm(cartesian(indAtom1, 2:end) - cartesian(indAtom2, 2:end));
z1 = cartesian(indAtom1, 1);
z2 = cartesian(indAtom2, 1);

coeffs = [ ...
    -2.630575471115972
    1.155281331403750
    0.406101379591945
    -0.006491787214362
    0.083551824390146
    -0.016819474606420
    -0.033933135725531
    0.000967977472307
    -0.292248868653409
    0.411978769668091
    -0.201582865554906
    0.317088796767966
    0.017315986848049
    0.001901847913666
    2.544990894929442];

criterions(1) = Feature(z1, z2) * coeffs;
criterions(2) = criterions(1) ./ 1.05 .* ...
    ~~max(min(MaxNumBonds(z1), MaxNumBonds(z2))-1, 0);
criterions(3) = criterions(1) ./ 1.15 .* ...
    ~~max(min(MaxNumBonds(z1), MaxNumBonds(z2))-2, 0);
criterions(1) = 1.2 .* criterions(1);

numBonds = 3 - sum(distance>criterions);

end

function feat = Feature(z1, z2)
[row1, col1] = FromAtomicNumber(z1);
[row2, col2] = FromAtomicNumber(z2);

feat = [row1+row2 row1^2+row2^2 row1*row2 row1^2*row2^2 ...
    col1+col2 col1^2+col2^2 col1*col2 col1^2*col2^2, ...
    row1*z1+row2*z2, sqrt(z1)+sqrt(z2) sqrt(z1*z2) z1+z2 z1^2+z2^2 z1*z2, ...
    1];
end

function [row, col, rowMax] = FromAtomicNumber(atomicNum)
if(atomicNum<=2)
    row = 1;
    col = atomicNum;
    rowMax = 2;
elseif(atomicNum<=10)
    row = 2;
    col = atomicNum-2;
    rowMax = 8;
elseif(atomicNum<=18)
    row = 3;
    col = atomicNum-10;
    rowMax = 8;
else
    row = 4;
    col = atomicNum-18;
    rowMax = 18;
end
end
