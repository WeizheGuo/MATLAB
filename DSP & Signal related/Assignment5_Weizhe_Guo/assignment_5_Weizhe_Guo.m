% Name: Weizhe Guo
% Date: 02/21/2018
% Assignment 5

load fisheriris
Flowerinst = cell(150,1);
for i = 1:150
    Flowerinst{i,1} = Flower(meas(i,1),meas(i,2),meas(i,3),meas(i,4),species{i,1});
end

n = randi(150);
testflower = Flower(meas(n,1),meas(n,2),meas(n,3),meas(n,4),species{n,1});
testlength = testflower.getSLength;
if testlength == meas(n,4)
    testresult = 1;
else
    testresult = 0;
end

testflower.report();