clear;
clc;
format long;
%Pi using a for loop
N = 1000;
count = 0;
for j = 1:N
    
    if (rand(1)^2 + rand(1)^2) < 1.0
        count = count + 1;
    end
    Pi(j) = 4*count/j;
end




plot([1:N],pi*ones(N,1));
hold on;
plot([1:N],Pi,'o');