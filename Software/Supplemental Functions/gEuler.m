%Input: complex: is a structure of finite dimensional simplicial complex (Here is for 3D images).
%       complex.V is n1x3 vertices
%       complex.E is n2x2 edges
%       complex.F is n3x3 faces
%       complex.T=[]; for 2D case.
%       fun: The function value on the vertices
%       stepsize: discrete Euler Characteristic curve into such dimensional
%       vector. (e.g. 100);
%Output: Stepsize_by_2 matrix. Euler Characteristics (second column) for scales of function values (first column).

function [kai]=gEuler(complex,fun,stepsize)
V=complex.V;%vertices:n1x3
E=complex.E;%edge:n2x2
F=complex.F;%face:n3x3
T=complex.T;%tetrahedra:n4x4
if length(V)~=length(fun)
    fprintf('The size of function should be same as the number of vertices');
    return
end
fe=zeros(size(E,1),1);
ff=zeros(size(F,1),1);
ft=zeros(size(T,1),1);
fe=max(fun(E)')';
ff=max(fun(F)')';
ft=max(fun(T)')';
%stepsize=10;
threshold=min(fun):(max(fun)-min(fun))/stepsize:max(fun);
kai=zeros(stepsize+1,2);
kai(:,1)=threshold';
for i=1:length(threshold)
    v=length(find(fun<=threshold(i)));
    e=length(find(fe<=threshold(i)));
    f=length(find(ff<=threshold(i)));
    t=length(find(ft<=threshold(i)));
    kai(i,2)=v-e+f-t;
end
