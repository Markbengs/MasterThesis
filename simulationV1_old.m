function r = SimulationV1_old(V,L,N,S,R0,RT,VD,tSpiral,tCircle,wc,draw) %args
%V = 13; R0 = 25; RT = 25; L = 30; N = 10; VD = 1/5; tStraight = 5; tSpiral = 30; tCircle = 30; S = 1; wc = [0;0;0];
mDrouge = 0.5; g = 9.82; e1 = [1;0;0]; e2 = [0;1;0]; e3 = [0;0;1];  rho = 1.225; d = 0.003; %d = 0.00635;

%Coefficients for drogue with 30cm diameter
CdD = 0.42; %0.42
CdL = 0.01; %0.01

L0 = L/N; mSegment = 0.007*L0; %0.0238
m = [mSegment*ones(1,N-1),mDrouge];
M = diag(repelem(m,3));
Fg = m*g.*e3;

h = 0.01;
tStraight = 8;

p = zeros(3,N+1); Dp = zeros(3,N+1);
for i=1:N+1
    p(:,i) = [-tStraight*V; R0; (i-1)*L0]; %[p1x dp1x p1y dp1y...dp1z p2x dp2x...
    Dp(:,i) = [V; 0; 0];
end

A = zeros(N,3*N);
b = zeros(N,1);

%gamma1 = 0.9; gamma2 = 0.036; seems to work
gamma1 = 0.9; gamma2 =  0.036; %gamma1 = 0.05; gamma2 = 0.002;

for q = 0:h:tStraight
    p(:,1) = [(q-tStraight)*V; R0; 0];
    Dp(:,1) = [V,0,0];
    
    A(1,1:3) = transpose(p(:,2)-p(:,1));
    b(1) = -norm(Dp(:,2)-Dp(:,1))^2;%+transpose(p(:,2)-p(:,1))*DDpm(:,q);
    for i = 2:N
        A(i,3*(i-2)+(1:3)) = -transpose(p(:,i+1)-p(:,i)); %l(:,i)
        A(i,3*(i-1)+(1:3)) = transpose(p(:,i+1)-p(:,i)); %l(:,i)
        b(i) = -norm(Dp(:,i+1)-Dp(:,i))^2;
    end
    
    lVector = p(:,2:N+1)-p(:,1:N);
    vVector = max(((Dp(:,1:N)-wc)+(Dp(:,2:N+1)-wc))/2, 0.00001);
    
    lScalar = sqrt(sum(lVector.^2));
    vScalar = sqrt(sum(vVector.^2));
    
    alpha = acos(sum(lVector.*vVector)./(lScalar.*vScalar));
    
    eD = -vVector./vScalar;
    eL = -cross(cross(vVector,lVector),vVector)./sqrt(sum(cross(cross(vVector,lVector),vVector).^2));
    
    Vp = vScalar.*cos(alpha);
    Vn = vScalar.*sin(alpha);
    
    Mp = Vp/343;
    Mn = Vn/343;
    
    Cf = 0.038-0.0425*Mp;
    Cn = 1.17+Mn/40-(Mn.^2)/4+5*(Mn.^3)/8;
    
    CD = Cf+Cn.*sin(alpha).^3;
    CL = Cn.*sin(alpha).^2.*cos(alpha);
    
    FsDrag = (rho*CD.*lScalar*d.*vScalar.^2.*eD)/2;
    FsLift = (rho*CL.*lScalar*d.*vScalar.^2.*eL)/2;
    
    edL = -cross(cross(Dp(:,end),e3),Dp(:,end))./sqrt(sum(cross(cross(Dp(:,end),e3),Dp(:,end))).^2);
    Fdrag = [(FsDrag(:,2:N)+FsDrag(:,1:N-1))/2,FsDrag(:,end)/2 - (rho*CdD*S*norm(Dp(:,end))*Dp(:,end))/2];
    Flift = [(FsLift(:,2:N)+FsLift(:,1:N-1))/2,(rho*CdL*S*norm(Dp(:,end))^2.*edL)/2 + FsLift(:,end)/2];
    
    F  = Fg+Fdrag+Flift;
    a = F./m;
    aVec = reshape(a,[3*N,1]);
    
    phi = lScalar.^2-L0^2;
    psi = sum(lVector.*(Dp(:,2:N+1)-Dp(:,1:N)));
    
    Dphi = zeros(3*N,N); %actually transpose of derivative
    Dpsi = zeros(3*N,N);
    
    for i = 1:N
        Dphi(i*3-2:i*3,i) = -2*lVector(:,i);
        Dphi(3*i+1:3*i+3,i) = 2*lVector(:,i);
        Dpsi(i*3-2:i*3,i) = Dp(:,i+1)-Dp(:,i);
        Dpsi(3*i+1:3*i+3,i) = -Dp(:,i+1)+Dp(:,i);
    end
    Dphi = Dphi(4:end,:);
    Dpsi = Dpsi(4:end,:);
    
    xddot = aVec+M^(-1/2)*lsqminnorm(double(A*M^(-1/2)),double(b-A*aVec))-gamma1*Dphi*transpose(phi)-gamma2*Dpsi*transpose(psi);
    %xddot = aVec+M^(-1/2)*pinv(A*M^(-1/2))*(b-A*aVec)-gamma1*Dphi*transpose(phi)+gamma2*Dpsi*transpose(psi);
    
    Dp = Dp + [zeros(3,1),h*reshape(xddot,[3,N])];
    p =  p + h*Dp;
    
    if draw == 1 && mod(q,0.01)==0
        plot3(p(1,:),p(2,:),p(3,:))
        axis([-tStraight*V 30 -R0-1 R0+1 -5 35])
        view(180,0)
        drawnow
    end
end

for q = 0:h:tSpiral
    p(:,1) = [(R0-(R0-RT)/tSpiral*q)*sin(V*q/(R0-(R0-RT)/tSpiral*q));(R0-(R0-RT)/tSpiral*q)*cos(V*q/(R0-(R0-RT)/tSpiral*q));0];
    Dp(:,1) = [(R0-RT)*sin(V*q*tSpiral/(q*R0-R0*tSpiral-q*RT))/tSpiral+R0*tSpiral*V*cos(V*q*tSpiral/(R0*tSpiral-q*R0+q*RT))/(R0*tSpiral-q*R0+q*RT);...
        (RT-R0)*cos(V*q*tSpiral/(R0*tSpiral-q*R0+q*RT))/tSpiral+R0*tSpiral*V*sin(V*q*tSpiral/(q*R0-q*RT-R0*tSpiral))/(R0*tSpiral-q*R0+q*RT); 0];
    DDpm = [R0^2*tSpiral^3*V^2*sin(V*q*tSpiral/(q*R0-q*RT-R0*tSpiral))/(R0*tSpiral-q*R0+q*RT)^3;...
        -R0^2*tSpiral^3*V^2*cos(V*q*tSpiral/(R0*tSpiral-q*R0+q*RT))/(R0*tSpiral-q*R0+q*RT)^3;0];
    
    A(1,1:3) = transpose(p(:,2)-p(:,1));
    b(1) = -norm(Dp(:,2)-Dp(:,1))^2+transpose(p(:,2)-p(:,1))*DDpm;
    for i = 2:N
        A(i,3*(i-2)+(1:3)) = -transpose(p(:,i+1)-p(:,i)); %l(:,i)
        A(i,3*(i-1)+(1:3)) = transpose(p(:,i+1)-p(:,i)); %l(:,i)
        b(i) = -norm(Dp(:,i+1)-Dp(:,i))^2;
    end
    
    lVector = p(:,2:N+1)-p(:,1:N);
    vVector = ((Dp(:,1:N)-wc)+(Dp(:,2:N+1)-wc))/2;
    
    lScalar = sqrt(sum(lVector.^2));
    vScalar = sqrt(sum(vVector.^2));
    
    alpha = acos(sum(lVector.*vVector)./(lScalar.*vScalar));
    
    eD = -vVector./vScalar;
    eL = -cross(cross(vVector,lVector),vVector)./sqrt(sum(cross(cross(vVector,lVector),vVector).^2));
    
    Vp = vScalar.*cos(alpha);
    Vn = vScalar.*sin(alpha);
    
    Mp = Vp/343;
    Mn = Vn/343;
    
    Cf = 0.038-0.0425*Mp;
    Cn = 1.17+Mn/40-(Mn.^2)/4+5*(Mn.^3)/8;
    
    CD = Cf+Cn.*sin(alpha).^3;
    CL = Cn.*sin(alpha).^2.*cos(alpha);
    
    FsDrag = (rho*CD.*lScalar*d.*vScalar.^2.*eD)/2;
    FsLift = (rho*CL.*lScalar*d.*vScalar.^2.*eL)/2;
    
    edL = -cross(cross(Dp(:,end),e3),Dp(:,end))./sqrt(sum(cross(cross(Dp(:,end),e3),Dp(:,end))).^2);
    Fdrag = [(FsDrag(:,2:N)+FsDrag(:,1:N-1))/2,FsDrag(:,end)/2 - (rho*CdD*S*norm(Dp(:,end))*Dp(:,end))/2];
    Flift = [(FsLift(:,2:N)+FsLift(:,1:N-1))/2,(rho*CdL*S*norm(Dp(:,end))^2.*edL)/2 + FsLift(:,end)/2];
    
    F  = Fg+Fdrag+Flift;
    a = F./m;
    aVec = reshape(a,[3*N,1]);
    
    phi = lScalar.^2-L0^2;
    psi = sum(lVector.*(Dp(:,2:end)-Dp(:,1:end-1)));
    
    Dphi = zeros(3*N,N); %actually transpose of derivative
    Dpsi = zeros(3*N,N);
    
    for i = 1:N
        Dphi(i*3-2:i*3,i) = -2*lVector(:,i);
        Dphi(3*i+1:3*i+3,i) = 2*lVector(:,i);
        Dpsi(i*3-2:i*3,i) = Dp(:,i+1)-Dp(:,i);
        Dpsi(3*i+1:3*i+3,i) = -Dp(:,i+1)+Dp(:,i);
    end
    Dphi = Dphi(4:end,:);
    Dpsi = Dpsi(4:end,:);
    
    xddot = aVec+M^(-1/2)*lsqminnorm(double(A*M^(-1/2)),double(b-A*aVec))-gamma1*Dphi*transpose(phi)-gamma2*Dpsi*transpose(psi);
    %xddot = aVec+M^(-1/2)*pinv(A*M^(-1/2))*(b-A*aVec)-gamma1*Dphi*transpose(phi)+gamma2*Dpsi*transpose(psi);
    
    Dp = Dp + [zeros(3,1),h*reshape(xddot,[3,N])];
    p =  p + h*Dp;
    
    if draw == 1 && mod(q,0.1)==0
        plot3(p(1,:),p(2,:),p(3,:))
        axis([-tStraight*V 30 -R0-1 R0+1 -5 35])
        %view(0,90)
        %view(90,0)
        drawnow
    end
end

tLastLap = round(2*pi*R0/V,2);
test = zeros(3,round(tLastLap/h));
for q = tSpiral:h:tCircle+tSpiral
    p(:,1) = [RT*sin(V*q/RT);RT*cos(V*q/RT); (q-tSpiral)*VD];
    Dp(:,1) = [V*cos(V*q/RT); -V*sin(V*q/RT); VD];
    DDpm = [-V^2*sin(V*q/RT)/RT; -V^2*cos(V*q/RT)/RT; 0];
    
    A(1,1:3) = transpose(p(:,2)-p(:,1));
    b(1) = -norm(Dp(:,2)-Dp(:,1))^2+transpose(p(:,2)-p(:,1))*DDpm;
    for i = 2:N
        A(i,3*(i-2)+(1:3)) = -transpose(p(:,i+1)-p(:,i)); %l(:,i)
        A(i,3*(i-1)+(1:3)) = transpose(p(:,i+1)-p(:,i)); %l(:,i)
        b(i) = -norm(Dp(:,i+1)-Dp(:,i))^2;
    end
    
    lVector = p(:,2:N+1)-p(:,1:N);
    vVector = ((Dp(:,1:N)-wc)+(Dp(:,2:N+1)-wc))/2;
    
    lScalar = sqrt(sum(lVector.^2));
    vScalar = sqrt(sum(vVector.^2));
    
    alpha = acos(sum(lVector.*vVector)./(lScalar.*vScalar));
    
    eD = -vVector./vScalar;
    eL = -cross(cross(vVector,lVector),vVector)./sqrt(sum(cross(cross(vVector,lVector),vVector).^2));
    
    Vp = vScalar.*cos(alpha);
    Vn = vScalar.*sin(alpha);
    
    Mp = Vp/343;
    Mn = Vn/343;
    
    Cf = 0.038-0.0425*Mp;
    Cn = 1.17+Mn/40-(Mn.^2)/4+5*(Mn.^3)/8;
    
    CD = Cf+Cn.*sin(alpha).^3;
    CL = Cn.*sin(alpha).^2.*cos(alpha);
    
    FsDrag = (rho*CD.*lScalar*d.*vScalar.^2.*eD)/2;
    FsLift = (rho*CL.*lScalar*d.*vScalar.^2.*eL)/2;
    
    edL = -cross(cross(Dp(:,end),e3),Dp(:,end))./sqrt(sum(cross(cross(Dp(:,end),e3),Dp(:,end))).^2);
    Fdrag = [(FsDrag(:,2:N)+FsDrag(:,1:N-1))/2,FsDrag(:,end)/2 - (rho*CdD*S*norm(Dp(:,end))*Dp(:,end))/2];
    Flift = [(FsLift(:,2:N)+FsLift(:,1:N-1))/2,(rho*CdL*S*norm(Dp(:,end))^2.*edL)/2 + FsLift(:,end)/2];
    
    F  = Fg+Fdrag+Flift;
    a = F./m;
    aVec = reshape(a,[3*N,1]);
    
    phi = lScalar.^2-L0^2;
    psi = sum(lVector.*(Dp(:,2:end)-Dp(:,1:end-1)));
    
    Dphi = zeros(3*N,N); %actually transpose of derivative
    Dpsi = zeros(3*N,N);
    
    for i = 1:N
        Dphi(i*3-2:i*3,i) = -2*lVector(:,i);
        Dphi(3*i+1:3*i+3,i) = 2*lVector(:,i);
        Dpsi(i*3-2:i*3,i) = Dp(:,i+1)-Dp(:,i);
        Dpsi(3*i+1:3*i+3,i) = -Dp(:,i+1)+Dp(:,i);
    end
    Dphi = Dphi(4:end,:);
    Dpsi = Dpsi(4:end,:);
    
    xddot = aVec+M^(-1/2)*lsqminnorm(double(A*M^(-1/2)),double(b-A*aVec))-gamma1*Dphi*transpose(phi)-gamma2*Dpsi*transpose(psi);
    %xddot = aVec+M^(-1/2)*pinv(A*M^(-1/2))*(b-A*aVec)-gamma1*Dphi*transpose(phi)+gamma2*Dpsi*transpose(psi);
    
    Dp = Dp + [zeros(3,1),h*reshape(xddot,[3,N])];
    p =  p + h*Dp;
    
    if draw == 1 && mod(q,0.1) == 0
        plot3(p(1,:),p(2,:),p(3,:))
        axis([-tStraight*V 30 -R0-1 R0+1 -5 35])
        %view(0,90)
        drawnow
    end
    
    if q > tCircle+tSpiral-tLastLap+1e-10
    %test(:,round((q-tSpiral)/h)+1) = p(:,end);
    test(:,round((q-tCircle-tSpiral+tLastLap)/h)) = p(:,end);
    end
end
plot3(test(1,:),test(2,:),test(3,:))
displacement = mean(test(1:2,:),2);
radius = sqrt(sum((displacement-test(1:2,:)).^2));
%r = mean(radius);
%r = mean(displacement(1:2));
r = displacement;
end