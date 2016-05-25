X=[1 2 3; 3 4 5];

%S=[X;X];

tp=[6 7 8];

theta=[1 0.5 0.25];

v1 = exp(-sum(theta.*(X(1,:) - tp).^2))
v2 = exp(-sum(theta.*(X(2,:) - tp).^2))
v3 = exp(-sum(theta.*(X(3,:) - tp).^2))

v1_1 = -2*theta(1)*(X(1,1)-tp(1,1))*v1
%v1_2 = -2*theta(2)*(X(1,2)-tp(1,2))*v1
%v1_3 = -2*theta(3)*(X(1,3)-tp(1,3))*v1

v2_1 = -2*theta(1)*(X(2,1)-tp(1,1))*v2
%v2_2 = -2*theta(2)*(X(2,2)-tp(1,2))*v2
%v2_3 = -2*theta(3)*(X(2,3)-tp(1,3))*v2

v3_1 = -2*theta(1)*(X(3,1)-tp(1,1))*v3
%v3_2 = -2*theta(2)*(X(3,2)-tp(1,2))*v3
%v3_3 = -2*theta(3)*(X(3,3)-tp(1,3))*v3