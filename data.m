%loading
load('X.dat')
load('Y.dat')
load('Z.dat')  
A=X( : ,1);
B=Y( : ,1);
C=Z( : ,1);
D=[ 1,0 ];
E=[ 1,1 ];
X = kron(X,E);
Y = kron(Y,E);
Z = kron(Z,D);
%preparing
X=[X A];
Y=[Y B];
Z=[Z C];