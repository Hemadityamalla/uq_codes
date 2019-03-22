clear;clc;format long;

for i=1:5 %functions
   for j=1:5 %quadrature
       fprintf("Function no. %i and Quadrature %i \n", i,j);
       Generalized_Implicit_Quad(j,i);
   end
end