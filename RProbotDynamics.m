clear all;
close all;
syms q1 q2 dq1 dq2 ddq1 ddq2 g m1 m2 I1 I2 L real

%STEP 1: Find O0,O1, Oc1 Oc2
O0 = [0; 0; 0];
O1 = [q2*cos(q1); q2*sin(q1); 0];
O2 = [(q2+2*L)*cos(q1); (q2+2*L)*sin(q1); 0];

Oc1 = [0; 0; 0];
Oc2 = [(q2+L)*cos(q1); (q2+L)*sin(q1); 0];

%STEP 2: Find Z0,Z1
Z0 = [0; 0; 1]; %third column of T00 because the rotation is along z axis
Z1 = [cos(q1); sin(q1); 0];%first collumn of T01 because the motion is along x axis

zeros = zeros(3,1);
%STEP 3: Find Jacobians for v1,v2,w1,w2
Jv1 = [cross(Z0,Oc1-O0), zeros];
Jv2 = [cross(Z0,Oc2-O0),Z1];

Jw1 = [Z0, zeros];
Jw2 = [Z0, zeros];

%rotation and translation matrices
R1 = [cos(q1) -sin(q1) 0
    sin(q1) cos(q1) 0
    0 0 1];

R2 = R1;

%STEP 4: Find kinetic energy
%generalized inertia matrix
B1 = simplify(m1*(Jv1')*Jv1 + Jw1'*R1'*I1*R1*Jw1); 
B2 = simplify(m2*(Jv2')*Jv2 + Jw2'*R2'*I2*R2*Jw2);
B = simplify(B1+B2);

%STEP 5: Find Coriolis ans centrifugal forces
q = [q1;q2];
dq = [dq1;dq2];
ddq = [ddq1;ddq2];
C = Coriolis(B,q,dq,2);
C = simplify(C);

%STEP 6: Find potential energy
U1 = 0;
U2 = m2*g*(q2+L)*sin(q1);
U = U1+U2;
%find the gravity
G1 = diff(U,q1);
G2 = diff(U,q2);
G = [G1; G2];

%STEP 7: Equation of motion
tau = B*ddq+C*dq+G;
ddq_num = [0.1; 0.1];
dq_num = [1; 1];
B_num = subs(B,{m1, m2, I1, I2, L, q1, q2},{2 2 1 2 0.2 0 1});
C_num = subs(C ,{m1, m2, I1, I2, L, q1, q2, dq1, dq2},{2 2 1 2 0.2 0 1 1 1});
G_num = subs(G ,{m1, m2, I1, I2, L, g, q1, q2},{2 2 1 2 0.2 9.81 0 1});
tau_num = B_num*ddq_num+C_num*dq_num+G_num;


%direct dynamics
B(q1,q2) = subs(B,{m1, m2, I1, I2, L},{2 2 1 2 0.2});
C(q1,q2,dq1,dq2) = subs(C*dq ,{m1, m2, I1, I2, L},{2 2 1 2 0.2});
G(q1,q2) = subs(G ,{m1, m2, I1, I2, L, g},{2 2 1 2 0.2 9.81});

q1_0 = pi/2;
q2_0 = 5;
dq1_0 = 0;
dq2_0 = 0;
dt = 0.1;

n = 50;
u = [0 0];

for i=1:n
%     u = [sin(2*i*dt); cos(i*dt)];
    u1p(i) = u(1);
    u2p(i) = u(2);
    q1p(i) = q1_0;
    q2p(i) = q2_0;
    dq1p(i) = dq1_0;
    dq2p(i) = dq2_0;
    ddq = inv(B(q1_0,q2_0))*(u-C(q1_0,q2_0,dq1_0,dq2_0)-G(q1_0,q2_0));
    
    %joint accelerations
    ddq1p(i) = ddq(1);
    ddq2p(i) = ddq(2);
    
    %joint velocities
    dq1_0 = dq1p(i) + double(ddq(1)*dt);
    dq2_0 = dq2p(i) + double(ddq(2)*dt);
    
    %joint position
    q1_0 = q1p(i) + dq1_0*dt;
    q2_0 = q2p(i) + dq2_0*dt;
    
end

t = 0:dt:dt*(n-1);

%position plot
figure
plot(t,q1p)
xlabel('time, sec')
ylabel('position, rad')
title('Position q1 graph')

figure
plot(t,q2p)
xlabel('time, sec');
ylabel('position, m')
title('Position q2 graph')

%velocity plot
figure
plot(t,dq1p)
xlabel('time, sec');
ylabel('velocity, rad/s')
title('Velocity q1 graph')

figure
plot(t,dq2p)
xlabel('time, sec')
ylabel('velocity, m/s')
title('Velocity q2 graph')

%acceleration plot
figure
plot(t,ddq1p)
xlabel('time, sec')
ylabel('acceleration, rad/s^2')
title('Acceleration q1 graph')

figure
plot(t,ddq2p)
xlabel('time, sec')
ylabel('acceleration, m/s^2')
title('Acceleration q2 graph')

figure
plot(t,u1p,t,u2p);
xlabel('time, sec');
title('Forces and torques graph');
legend('u_1','u_2')

function C = Coriolis(B,q,dq,n)
syms C
for k = 1:n
    for j = 1:n
        C(k,j) = sym(0);
        for i = 1:n
            c_ijk = 0.5*(diff(B(k,j),q(i))+diff(B(k,i),q(j))-diff(B(i,j),q(k)));
            C(k,j) = C(k,j) + c_ijk*dq(i);
        end
    end
end
end
