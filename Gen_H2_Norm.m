function mu = Gen_H2_Norm(A,B,C)
% MAE 598 - Natalie Beaulieu - Dr. Matthew Peet
% Function solves LMI for Generalized H2 Norm
% A must be Hurwitz
G = tf(A,B,C,0);
n = size(A,1);
rho = 0.001;
options = sdpsettings('verbose',0,'solver','sedumi');

% LMI one
% Variables
P1 = sdpvar(n);
mu1 = sdpvar(n);

% Constraints
F1 = [P1 >= rho*eye(n)];
F1 = [F1, (A'*P1)+(P1*A) <= -rho*eye(n)];
F1 = [F1, (P1*B) <= -rho*eye(n)];
F1 = [F1, -mu1 <= -rho*eye(n)];
F1 = [F1, C' >= rho*eye(n)];

% Optimize
optimize(F1,[],options);
mu1 = value(mu1);

% LMI two
% Variables
Q2 = sdpvar(n);
mu2 = sdpvar(n);

% Constraints
F2 = [Q2 >= rho*eye(n)];
F2 = [F2, (Q2*A')+(A*Q2) <= -rho*eye(n)];
F2 = [F2, (B) <= -rho*eye(n)];
F2 = [F2, -mu2 <= -rho*eye(n)];
F2 = [F2, (Q2*C') >= rho*eye(n)];

% Optimize
optimize(F2,[],options);
mu2 = value(mu2);

% LMI three
% Variables
P3 = sdpvar(n);
V3 = sdpvar(n);
mu3 = sdpvar(n);

% Constraints
F3 = [P3 >= rho*eye(n)];
F3 = [F3, -(V+V') <= -rho*eye(n)];
F3 = [F3, ((V'*A)+P) <= -rho*eye(n)];
F3 = [F3, (V'*B) <= -rho*eye(n)];
F3 = [F3, V' <= -rho*eye(n)];
F3 = [F3, -mu3 <= -rho*eye(n)];
F3 = [F3, (C') >= rho*eye(n)];

% Optimize
optimize(F3,[],options);
mu3 = value(mu3);

mu = [mu1,mu2,mu3];
end

