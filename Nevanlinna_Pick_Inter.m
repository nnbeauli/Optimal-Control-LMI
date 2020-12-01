function [gamma] = Nevanlinna_Pick_Inter(lambda, U, V)
% MAE 598 - Natalie Beaulieu - Dr. Matthew Peet
% Function solves LMI for Nevanlinna Pick Interpolation with scaling using bisection
% Lambda is sequence of system in complex plane,
% U and V are vectors of target points,
% D is set of m x m block-diagonal matrices
% A is diagonal matrix of lambdas

%System
A = diag(lambda);
rho = 0.0001;
options = sdpsettings('verbose',0,'solver','mosek');

% Variables
Gin = sdpvar(size(A,1),size(A,2));
Gout=sdpvar(size(A,1),size(A,2));
D = sdpvar(size(lambda,1),size(lambda,2));
P = sdpvar(size(U,1));
gamma =sdpvar(1);
gmmi = gamma^2;

% Matrices
mat1 = [(Im(A)*Gin)+(Gin*A)-(Im(U)*P*U)];
mat2 = [(Im(A)*Gout)+(Gout*A)-(Im(V)*P*V)];
mat3 = [(gmmi*Gin)-Gout];

% Constraints
F = [D == Im(D)];
F = [F, D >= rho*eye(size(D,1),size(D,2))];
F = [F, P == Im(D)*D];
F = [F, P >= rho*eye(size(P,1),size(P,2))];
F = [F, mat1 == rho*eye(size(mat1,1),size(mat1,2))];
F = [F, mat2 == rho*eye(size(mat2,1),size(mat2,2))];
F = [F, mat3 >= rho*eye(size(mat3,1),size(mat3,2))];

% Optimization
bisection(F,gmmi,options)
gmmi = value(gmmi);
gamma = sqrt(gmmi);
disp('Optimal gamma given by Navenlinna Pick Interpolation is:')
disp(gmmi)
end