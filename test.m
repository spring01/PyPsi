% A simple test script creates a MatPsi object named matpsi 

cart = [ ...
    8 0.0 0.0 0.0;
    1 0.0 0.0 1.0;
    1 0.0 1.0 0.0;];
matpsi = MatPsi2(cart, '6-31g*');

