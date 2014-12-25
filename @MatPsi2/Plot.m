function Plot(obj)
cartesian = obj.Molecule_Geometry() .* 0.529177249;

numAtoms = size(cartesian, 1);
maxNumBondsVec = zeros(numAtoms, 1);
for iAtom = 1:numAtoms
    maxNumBondsVec(iAtom) = MaxNumBonds(cartesian(iAtom, 1));
end

[maxNumBondsVec, order] = sort(maxNumBondsVec, 'descend');
sorted_cartesian = cartesian(order,:);

numBondsMat = zeros(numAtoms);
for iAtom = 1:numAtoms
    for jAtom = iAtom+1:numAtoms
        numBondsMat(iAtom, jAtom) = NumberOfBonds(sorted_cartesian, iAtom, jAtom);
    end
end
numBondsMat = numBondsMat + numBondsMat';

for iAtom = 1:numAtoms
    [temp, tempOrder] = sort(numBondsMat(iAtom, :));
    cumTemp = cumsum(temp);
    for i = 1:length(cumTemp)
        if(cumTemp(i) > maxNumBondsVec(iAtom))
            numBondsMat(iAtom, tempOrder(i)) = maxNumBondsVec(iAtom) - cumTemp(i-1);
            numBondsMat(tempOrder(i), iAtom) = numBondsMat(iAtom, tempOrder(i));
            numBondsMat(iAtom, tempOrder(i+1:end)) = 0;
            numBondsMat(tempOrder(i+1:end), iAtom) = 0;
            break;
        end
    end
end

bondLines = {};
for iAtom = 1:numAtoms
    iPos = sorted_cartesian(iAtom, 2:end);
    for jAtom = iAtom+1:numAtoms
        jPos = sorted_cartesian(jAtom, 2:end);
        numBonds = numBondsMat(iAtom, jAtom);
        vector = iPos - jPos;
        normal_vector = ones(1, 3);
        normal_vector(3) = -(vector(1) + vector(2))/vector(3);
        normal_vector = normal_vector ./ norm(normal_vector);
        switch(numBonds)
            case(1) % single bond
                bondLines{end+1} = [iPos; jPos];
            case(2) % double bond
                bondLines{end+1} = [iPos; jPos] ...
                    - 0.03.*[normal_vector; normal_vector];
                bondLines{end+1} = [iPos; jPos] ...
                    + 0.03.*[normal_vector; normal_vector];
            case(3)% triple bond
                bondLines{end+1} = [iPos; jPos];
                bondLines{end+1} = [iPos; jPos] ...
                    - 0.05.*[normal_vector; normal_vector];
                bondLines{end+1} = [iPos; jPos] ...
                    + 0.05.*[normal_vector; normal_vector];
        end
    end
end

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

shapes = allShapes(sorted_cartesian(:, 1));
colors = allColors(sorted_cartesian(:, 1), :);
scatter3(sorted_cartesian(:, 2), sorted_cartesian(:, 3), sorted_cartesian(:, 4), shapes, colors, 'fill');

end


function numBonds = NumberOfBonds(cartesian, indAtom1, indAtom2)

distance = norm(cartesian(indAtom1, 2:end) - cartesian(indAtom2, 2:end));
z1 = cartesian(indAtom1, 1);
z2 = cartesian(indAtom2, 1);

coeffs = [-2.6306 1.1553 0.4061 -0.0065 0.0836 -0.0168 -0.0339, ...
    0.0010 -0.2922 0.4120 -0.2016 0.3171 0.0173 0.0019 2.5450]';

criterions(1) = Feature(z1, z2) * coeffs;
criterions(2) = criterions(1) ./ 1.1 .* ...
    ~~max(min(MaxNumBonds(z1), MaxNumBonds(z2))-1, 0);
criterions(3) = criterions(1) ./ 1.2 .* ...
    ~~max(min(MaxNumBonds(z1), MaxNumBonds(z2))-2, 0);
criterions(1) = criterions(1).*1.2;

numBonds = 3 - sum(distance>criterions);

end

function maxNumBonds = MaxNumBonds(atomicNum)
[~, col, rowMax] = FromAtomicNumber(atomicNum);
maxNumBonds = min(col, rowMax-col);
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
