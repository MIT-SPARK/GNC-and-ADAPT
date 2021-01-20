function flag = areBinaryWeights(weights)
% Check the input vector is (approximately) a binary vectory
flag = true;
th = 1e-12;
for i = length(weights):-1:1
    if weights(i) > th && abs(1 - weights(i)) > th
        flag = false;
        break
    end
end
end