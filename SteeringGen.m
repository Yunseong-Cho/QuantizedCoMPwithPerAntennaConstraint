function [a] = SteeringGen(theta, NumAntenna, NumPath)

a = zeros(NumAntenna, NumPath);

base = 0:1:NumAntenna-1';

for i = 1:NumPath
    a(:,i) = (1/sqrt(NumAntenna))*exp(-1j*pi*theta(i)*base);
end

end

