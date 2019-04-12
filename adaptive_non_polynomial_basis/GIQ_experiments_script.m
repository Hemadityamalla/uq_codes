clear;clc;format long;

for i=1:6 %functions
   for j=1:6 %quadrature
       fprintf("Function no. %i and Quadrature %i \n", i,j);
       Generalized_Implicit_Quad(j,i);
   end
end