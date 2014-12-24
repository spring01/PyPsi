function Plot(obj)

cartesian = [obj.Molecule_AtomicNumbers()' obj.Molecule_Geometry() .* 0.529177249];

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

geom = cartesian(:, 2:end);
shapes = allShapes(cartesian(:,1));
colors = allColors(cartesian(:,1), :);
bondLines = {};

for iAtom = 1:size(cartesian,1)
    % determine bonds
    for jAtom = iAtom+1:size(cartesian,1)
        numBonds = NumberOfBonds(cartesian, iAtom, jAtom);
        vector = geom(iAtom, :) - geom(jAtom, :);
        normal_vector = ones(1, 3);
        normal_vector(3) = -(vector(1) + vector(2))/vector(3);
        normal_vector = normal_vector ./ norm(normal_vector);
        switch(numBonds)
            case(1) % single bond
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)];
            case(2) % double bond
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)] ...
                    - 0.03.*[normal_vector; normal_vector];
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)] ...
                    + 0.03.*[normal_vector; normal_vector];
            case(3)% triple bond
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)];
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)] ...
                    - 0.05.*[normal_vector; normal_vector];
                bondLines{end+1} = [geom(iAtom, :); geom(jAtom, :)] ...
                    + 0.05.*[normal_vector; normal_vector];
        end
    end
end
figure();
axis vis3d;
hold;
for oneLine = bondLines
    line(oneLine{1}(:, 1), oneLine{1}(:, 2), oneLine{1}(:, 3), 'LineWidth', 4);
end
scatter3(geom(:,1), geom(:,2), geom(:,3), shapes, colors, 'fill');

end


function numBonds = NumberOfBonds(cartesian, atom1, atom2)

distance = norm(cartesian(atom1, 2:end) - cartesian(atom2, 2:end));
z1 = cartesian(atom1, 1);
z2 = cartesian(atom2, 1);

coeffs = [-2.6306 1.1553 0.4061 -0.0065 0.0836 -0.0168 -0.0339, ...
    0.0010 -0.2922 0.4120 -0.2016 0.3171 0.0173 0.0019 2.5450]';

criterions(1) = Feature(z1, z2) * coeffs;
criterions(2) = criterions(1) ./ 1.1 .* DoubleBondPossible(z1, z2);
criterions(3) = criterions(1) ./ 1.2 .* TripleBondPossible(z1, z2);
criterions(1) = criterions(1) .* 1.2;

numBonds = 3 - sum(distance>criterions);

end

function possible = DoubleBondPossible(z1, z2)
list = [4 5 6 7 8 12 13 14 15 16];
possible = ~~(sum(z1==list).*sum(z2==list));
end

function possible = TripleBondPossible(z1, z2)
list = [5 6 7 13 14 15];
possible = ~~(sum(z1==list).*sum(z2==list));
end

function feat = Feature(z1, z2)
[row1, col1] = FromElement(z1);
[row2, col2] = FromElement(z2);

feat = [row1+row2 row1^2+row2^2 row1*row2 row1^2*row2^2 ...
    col1+col2 col1^2+col2^2 col1*col2 col1^2*col2^2, ...
    row1*z1+row2*z2, sqrt(z1)+sqrt(z2) sqrt(z1*z2) z1+z2 z1^2+z2^2 z1*z2, ...
    1];
end

function [row, col] = FromElement(element)
if(element<=2)
    row = 1;
    col = element;
elseif(element<=10)
    row = 2;
    col = element-2;
elseif(element<=18)
    row = 3;
    col = element-10;
else
    row = 4;
    col = element-18;
end
end
