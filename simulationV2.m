function [r, varargout] = SimulationV2(V,L,N,S,R,VD,tCircle,wc,draw) %args
%V = 13; L = 30; N = 10; S = 1; R0 = 25; RT = 25; VD = 1/5; tSpiral = 20; tCircle = 40; wc = [0;0;0];
rho = 1.225; d = 0.003;  g = 9.82; e3 = [0;0;1]; %k = 6.3338e+04; %d = 0.00635; k=A*E/L0

L0 = L/N; mDrouge = 0.5; %wc = [0;0;0];
mSegment = 0.007*L0; %0.0238
m = [mSegment*ones(1,N-1),mDrouge];
Fg = m*g.*e3;

tStraight=8;

p = zeros(3,N+1); Dp=zeros(3,N+1);
for i=1:N+1
    p(:,i) = [-tStraight*V; R; (i-1)*L0]; %[p1x dp1x p1y dp1y...dp1z p2x dp2x...
    Dp(:,i) = [V; 0; 0];
end

%Coefficients for drogue with 30cm diameter
CdD = 0.42; %0.42
CdL = 0.01; %0.01

E = 1.5e+09; %1-1.5 seems to be right for polyester, Diolen
A = d^2*pi/4;
    
h = 0.002;
ttest = zeros(3,tCircle/h);
for i=0:h:tStraight
    
    pm = [V*(i-tStraight); R; 0];
    p(:,1) = pm;
    
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
    
    Fdrag = (FsDrag(:,2:N)+FsDrag(:,1:N-1))/2;
    Flift = (FsLift(:,2:N)+FsLift(:,1:N-1))/2;
    
    edL = -cross(cross(Dp(:,end),e3),Dp(:,end))./sqrt(sum(cross(cross(Dp(:,end),e3),Dp(:,end))).^2);
    Fdrag = [Fdrag, FsDrag(:,end)/2 - (rho*CdD*S*sqrt(sum(Dp(:,end).^2))*Dp(:,end))/2];
    Flift = [Flift, (rho*CdL*S*sum(Dp(:,end).^2).*edL)/2 + FsLift(:,end)/2];
    
    %Newton law model
    eT = lVector./lScalar;
    T = -E*A*(lScalar-L0).*eT/L0;
    Fs = T-[T(:,2:N),zeros(3,1)];
    
    F = Fdrag+Flift+Fs+Fg;
    a = F./m;
    
    Dp = Dp + [zeros(3,1),h*a];
    p =  p + h*Dp;
    
    if draw == 1 && mod(i,0.1) == 0
        plot3(p(1,:),p(2,:),p(3,:))
        set(gca,'Ydir','reverse','Zdir','reverse')
        axis([-70 30 -30 30 -5 35])
        %view(180,0)
        drawnow
    end
    ttest(:,round(i/h)+1) = p(:,end);
end

% r = (R0-RT)/tSpiral;
% for i=0:h:tSpiral
%     p(:,1) = [(R0-r*i)*sin(V*i/(R0-r*i));(R0-r*i)*cos(V*i/(R0-r*i));0]; %Reversed since NED coord system
%     
%     lVector = p(:,2:N+1)-p(:,1:N);
%     vVector = ((Dp(:,1:N)-wc)+(Dp(:,2:N+1)-wc))/2;
%     
%     lScalar = sqrt(sum(lVector.^2));
%     vScalar = sqrt(sum(vVector.^2));
%     
%     alpha = acos(sum(lVector.*vVector)./(lScalar.*vScalar));
%     
%     eD = -vVector./vScalar;
%     eL = -cross(cross(vVector,lVector),vVector)./sqrt(sum(cross(cross(vVector,lVector),vVector).^2));
%     
%     Vp = vScalar.*cos(alpha);
%     Vn = vScalar.*sin(alpha);
%     
%     Mp = Vp/343;
%     Mn = Vn/343;
%     
%     Cf = 0.038-0.0425*Mp;
%     Cn = 1.17+Mn/40-(Mn.^2)/4+5*(Mn.^3)/8;
%     
%     CD = Cf+Cn.*sin(alpha).^3;
%     CL = Cn.*sin(alpha).^2.*cos(alpha);
%     
%     FsDrag = (rho*CD.*lScalar*d.*vScalar.^2.*eD)/2;
%     FsLift = (rho*CL.*lScalar*d.*vScalar.^2.*eL)/2;
%     
%     Fdrag = (FsDrag(:,2:N)+FsDrag(:,1:N-1))/2;
%     Flift = (FsLift(:,2:N)+FsLift(:,1:N-1))/2;
%     
%     edL = -cross(cross(Dp(:,end),e3),Dp(:,end))./sqrt(sum(cross(cross(Dp(:,end),e3),Dp(:,end))).^2);
%     Fdrag = [Fdrag, FsDrag(:,end)/2 - (rho*CdD*S*sqrt(sum(Dp(:,end).^2))*Dp(:,end))/2];
%     Flift = [Flift, (rho*CdL*S*sum(Dp(:,end).^2).*edL)/2 + FsLift(:,end)/2];
%     
%     %Newton law model
%     eT = lVector./lScalar;
%     E = 1.5e+09; %1-1.5 seems to be right for polyester, Diolen
%     A = d^2*pi/4;
%     T = -E*A*(lScalar-L0).*eT/L0;
%     Fs = T-[T(:,2:N),zeros(3,1)];
%     
%     F = Fdrag+Flift+Fs+Fg;
%     a = F./m;
%     
%     Dp = Dp+[zeros(3,1),h*a];
%     p =  p + h*Dp;
%     
%     if draw == 1 && mod(i,0.1) == 0
%         plot3(p(1,:),p(2,:),p(3,:))
%         axis([-30 30 -30 30 -5 35])
%         %view(180,0)
%         drawnow
%     end
%     ttest(:,round(i/h)+1) = p(:,end);
% end

%test = zeros(3,tCircle/h);
for i=0:h:tCircle
    
    p(:,1) = [R*sin(V*i/R);R*cos(V*i/R); i*VD];
    
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
    
    Fdrag = (FsDrag(:,2:N)+FsDrag(:,1:N-1))/2;
    Flift = (FsLift(:,2:N)+FsLift(:,1:N-1))/2;
    
    edL = -cross(cross(Dp(:,end),e3),Dp(:,end))./sqrt(sum(cross(cross(Dp(:,end),e3),Dp(:,end))).^2);
    Fdrag = [Fdrag, FsDrag(:,end)/2 - (rho*CdD*S*sqrt(sum(Dp(:,end).^2))*Dp(:,end))/2];
    Flift = [Flift, (rho*CdL*S*sum(Dp(:,end).^2).*edL)/2 + FsLift(:,end)/2];
    
    %Newton law model
    eT = lVector./lScalar;
    T = -E*A*(lScalar-L0).*eT/L0;
    Fs = T-[T(:,2:N),zeros(3,1)];
    
    F = Fdrag+Flift+Fs+Fg;
    a = F./m;
    
    Dp = Dp+[zeros(3,1),h*a];
    p =  p + h*Dp;
    
    if draw == 1 && mod(i,0.1) == 0
        plot3(p(1,:),p(2,:),p(3,:))
        set(gca,'Ydir','reverse','Zdir','reverse')
        axis([-30 30 -30 30 -5 35])
        %view(180,0)
        drawnow
    end
    %test(:,round((i-tSpiral)/h)+1) = p(:,end);
    ttest(:,round(i/h)+1) = p(:,end);
end

tLastLap = round(2*pi*R/V,2);
test = ttest(:,end-tLastLap/h:end);
%plot3(test(1,:),test(2,:),test(3,:))
plot3(ttest(1,:),ttest(2,:),ttest(3,:))
displacement = mean(test(1:2,:),2);
radius = sqrt(sum((displacement-test(1:2,:)).^2));
r = mean(radius);
varargout{1} = displacement;
zfix = test(3,:)-(0:h:tLastLap)*VD;
varargout{2} = max(zfix)-min(zfix);
end