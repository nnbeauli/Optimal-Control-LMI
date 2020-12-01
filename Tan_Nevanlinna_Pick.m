function N = Tan_Nevanlinna_Pick(lambda,U,v)
% MAE 598 - Natalie Beaulieu - Dr. Matthew Peet
% Function solves LMI for Tangential Nevanlinna Pick Interpolation
% Lambda is sequence of system in open disk,
% U and V are vectors of target points of size m,
% A is diagonal matrix of lambdas

%System
A = diag(lambda);
n = size(A,1);
rho = 0.0001;
options = sdpsettings('verbose',0,'solver','mosek');

% Variables
N = sdpvar(n);

% Constraints
F = [(Im(A)*N)+(N*A)-((Im(U)*U)-(Im(V)*V)) == rho*eye(n)];
N = value(N);

if N>= 0 
    disp('N>=0 and thus tangential Nevanlinna-Pick interpolation can be solved to find analytic function H')
else
    disp('N is not >=0 and thus tangential Nevanlinna-Pick interpolation can not be solved to find analytic function H')
end
end

