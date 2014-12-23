% A simple test script that creates a MatPsi object named matpsi 

str1 = ['O',char(10),'H 1 R', char(10),'H 1 R 2 A',char(10),char(10),'R = .9', char(10),'A = 104.5'];
str2 = '6-31g*';
matpsi = MatPsi(str1,str2);

