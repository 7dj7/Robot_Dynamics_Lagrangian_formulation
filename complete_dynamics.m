clear; clc;

%% homogeneous transforms

n = 7; % DOF

% DH parameters symbols
q = sym('q', [n 1], 'real'); % generalized coordinates (joint angles)
d = sym('d', [n 1], 'real'); % link offsets
syms a1
syms g

% initial conditions for the configuration of Sawyer shown in the Figure.

q0 = [0 3*pi/2 0 pi 0 pi 3*pi/2];
d0 = [317 192.5 400 168.5 400 136.3 133.75];
a10 = 81;

% cell array of your homogeneous transformations; each Ti{i} is a 4x4 symbolic transform matrix
% NOTE: for symbolic arrays: q(1) = q1, q(2) = q2, etc.
Ti = cell(n,1);

% Homogeneous transformations solution
DH = [a1 -pi/2 d(1) q(1); 0 -pi/2 d(2) q(2); 0 -pi/2 d(3) q(3); 0 -pi/2 d(4) q(4); 0 -pi/2 d(5) q(5); 0 -pi/2 d(6) q(6); 0 0 d(7) q(7)]; % DH Parameter Matrix
T = eye(4);
for i = 1:n
    temp = compute_dh_matrix(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    T = T*temp;
    Ti{i} = T;
end

%% angular velocity jacobian (Jw)
% Initialize angular velocity jacobian as an nx1 cell array where each element is
% an 3xn symbolic matrix
Jw = arrayfun(@(x) sym(['Jw' num2str(x)], [3,n], 'real'), 1:n, 'UniformOutput', 0)';
% Angular velocity jacobian solution
Jw{1} = [[0;0;1] repmat([0;0;0],1,n-1)];
for i=2:n
       jw = [[0;0;1]];
       for j=1:i-1
           jw = [jw Ti{j}(1:3,3)];
       end
       jw = [jw repmat([0;0;0],1,n-i)];
       Jw{i} =jw;
end

%% linear velocity jacobian (Jv)
% the center of mass of each link measured relative to the link fixed frame
% like Ti and Jw, c is an nx1 cell array where each element is a symoblic vector/matrix
% for example: c{3} = [c3x c3y c3z]' is the center of mass of link 3 measured relative to frame 3
c = arrayfun(@(x) [sym(['c' num2str(x) 'x'], 'real'), sym(['c' num2str(x) 'y'], 'real'), ...
    sym(['c' num2str(x) 'z'], 'real')]', 1:n, 'UniformOutput', 0)';
% as with the angular velocity jacobian, the linear velocity jacobian is comprised of n 3xn
% symbolic matrices stored in a cell array. Jv{i} is the 3xn angular velocity jacobian of link i
Jv = cell(n,1);
% Linear velocity jacobian solution
P = eye(4);
for i=1:n    
P = Ti{i}*[[1;0;0;0] [0;1;0;0] [0;0;1;0] [c{i};1]];
    x = P(1,4);
    y = P(2,4);
    z = P(3,4);
    for j=1:n
        Jv{i} = [Jv{i} [diff(x,q(j));diff(y,q(j));diff(z,q(j))]];
    end        
end

%% potential energy
m = sym('m', [n 1], 'real'); % mass of each link

% Potential energy solution
P = eye(4);
PE = 0;

%% inertial matrix and kinetic energy
qd = sym('qd', [n 1], 'real'); % "q dot" - the first derivative of the q's in time (joint velocities)

% inertia tensor for each link relative to the inertial frame stored in an nx1 cell array
I = arrayfun(@(x) inertia_tensor(x), 1:n, 'UniformOutput', 0)';

% D = Inertia matrix solution & P = Poterntial Energy
D=0;
for i = 1:n
    D = D + (m(i)*Jv{i}'*Jv{i} + Jw{i}'*I{i}*Jw{i});
    P = Ti{i}*[[1;0;0;0] [0;1;0;0] [0;0;1;0] [c{i};1]];
    PE = PE + m(i)*g*P(3,4);
end

% KE = Kinetic energy solution
KE = 0.5*qd'*D*qd;
%% equations of motion

qdd = sym('qdd', [n 1], 'real'); % "q double dot" - the second derivative of the q's in time (joint accelerations)

% The Christoffel symbols
c = zeros(n,n,n,'sym');
for k = 1:n
    for i = 1:n
        for j =1:n
            c(i,j,k) = 0.5 * (diff(D(k,j),q(i)) + diff(D(k,i),q(j)) - diff(D(i,j),q(k)));
        end
    end
end

% The coriolis matrix
C = zeros(n,n,'sym');
for k = 1:n
    for j = 1:n
        temp = 0;
        for i = 1:n
            temp = temp + c(i,j,k)*qd(i);
        end
        C(j,k) = temp;
    end
end

% The gravitation terms
G = zeros(n,1,'sym');
for k = 1:n
    G(k) = diff(PE,q(k));
end

% eom_lhs = equations of motion
eom_lhs = D*qdd + C*qd + G;


