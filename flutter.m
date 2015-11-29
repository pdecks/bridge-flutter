% Patricia Decker
% CE 291F Project

% Solution of 2 DOF Flutter Equations (nondimensionalized)
% From "State-of-the-Art Methods for Calculating Flutter, Vortex-Induced,
% and Buffeting Responses of Bridge Structures" FHWA-RD-80-050 (1981)

clear all; close all; clc;

% Assumed/Derived Parameters
tran.w= 2*pi*0.1; % rad/sec, frequency of deck, heaving
tor.w = 2*tran.w; % rad/sec, frequency of deck, pitching

tran.xi= 0.01; % mechanical damping (% critical)
tor.xi= 0.01; 

% density of air
param.rho = 1.225/(515.138); % (lbf*s^2/ft)/ft^2 (density, english)

param.B = 100; % feet, deck width
    %param.B = 30.5; % meters, deck width
param.L = 4000; % feet, span
    %param.L = 1220; % meters, span

% coefficients due to generalized coordinates (see p29), in this case, for
% half sine wave vertical and torsional modes
param.C11 = param.L/2;
param.C22 = param.L/2;
param.C12 = param.L/2;

param.M = 711.8; % lb*sec^2/ft^2, constant mass of the deck/unit span
    %param.M = 3481; % kg*sec^2/m^2, constant mass of the deck/unit span
param.M1 = param.C11*param.M;

param.I = 857000; % lb*sec^2, constant moment of inertia/unit span
    %param.I = 3.9*10^5; % kg*sec^2, constant moment of inertia/unit span
param.I1 = param.C22*param.I;

% % values to check from p39
% (rho*(B^4)/I1)*C11
% (rho*(B^2)/M1)*C11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of the Flutter Equations using experimental flutter derivatives%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flutter Derivatives from Table and Corresponding Reduced Frequency
param.npoints = 6;
fder.A(1).graph = [0 0 0.75 0.70 0.68 0.70];
fder.A(2).graph = [0 -0.03 -0.05 -0.10 -0.14 -0.16];
fder.A(3).graph = [0 0 0.50 1.00 1.46 1.69];
fder.H(1).graph = [-0.67 -1.5 -2.05 -3.25 -4.25 -5.5];
fder.H(2).graph = [0 0 0.7 2.25 4.25 8.90];
fder.H(3).graph = [0 -0.05 -1.25 -3.35 -4.00 -5];
fder.U_NB.graph = [2 4 6 8 10 12];

% add (0,0) point for plotting
fder.U_NB0.graph = [0 fder.U_NB.graph]; fder.U_NB0.name = 'U/NB'; fder.U_NB0.label = 'U/NB = 2*pi/K';
fder.A0(1).graph = [0 fder.A(1).graph]; fder.A0(1).name = 'A_1*'; fder.A0(1).label = 'A_1*';
fder.A0(2).graph = [0 fder.A(2).graph]; fder.A0(2).name = 'A_2*'; fder.A0(2).label = 'A_2*';
fder.A0(3).graph = [0 fder.A(3).graph]; fder.A0(3).name = 'A_3*'; fder.A0(3).label = 'A_3*';
fder.H0(1).graph = [0 fder.H(1).graph]; fder.H0(1).name = 'H_1*'; fder.H0(1).label = 'H_1*';
fder.H0(2).graph = [0 fder.H(2).graph]; fder.H0(2).name = 'H_2*'; fder.H0(1).label = 'H_2*';
fder.H0(3).graph = [0 fder.H(3).graph]; fder.H0(3).name = 'H_3*'; fder.H0(1).label = 'H_3*';
fder.xaxis.graph = [0 0 0 0 0 0 0];

% plot of aerodynamic coefficients/flutter derivatives from Scanlan 1981
figure
subplot(3,2,1);plot(fder.U_NB0.graph,fder.A0(1).graph,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('A_1*');axis([0 12 0 3]);grid on
subplot(3,2,2);plot(fder.U_NB0.graph,-fder.H0(1).graph,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('-H_1*');axis([0 12 0 6]);grid on
subplot(3,2,3);plot(fder.U_NB0.graph,fder.A0(2).graph,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('A_2*');axis([0 12 -0.2 0.4]);grid on
subplot(3,2,4);plot(fder.U_NB0.graph,fder.H0(2).graph,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('H_2*');axis([0 12 -2 8]);grid on
subplot(3,2,5);plot(fder.U_NB0.graph,fder.A0(3).graph,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('A_3*');axis([0 12 0 3]);grid on
subplot(3,2,6);plot(fder.U_NB0.graph,-fder.H0(1).graph,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('-H_3*');axis([0 12 0 8]);grid on

%% FIND ROOTS OF COMPLEX POLYNOMIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cubic equation
% b3*x^3 + b2*x^2 + b1x + b0 = 0
j = 1;
for j = 1:param.npoints
    cubic.b0(j) = 2*tran.xi*((tor.w/tran.w)^2) + 2*tor.xi*tor.w/tran.w;
    cubic.b1(j) = -(param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).graph(j)*((tor.w/tran.w)^2) - (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).graph(j);
    cubic.b2(j) = -2*tor.xi*tor.w/tran.w - 2*tran.xi - 2*tran.xi*(param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).graph(j);
    cubic.b3(j) = (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).graph(j) + (param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).graph(j) + (param.rho*(param.B^2)/param.M1)*(param.rho*(param.B^4)/param.I1)*(param.C11*param.C22*fder.H(1).graph(j)*fder.A(3).graph(j) - (param.C12^2)*fder.A(1).graph(j)*fder.H(3).graph(j));
    cubic.roots_b(:,1) = roots([cubic.b3(j) cubic.b2(j) cubic.b1(j) cubic.b0(j)]);
    k = 1;
    for k=1:3
        if cubic.roots_b(k,1)>0
            cubic.root_b(j) = cubic.roots_b(k,1);
        end % end of positive root if statement
    end % end of positive-root-finding loop
    % assemble points for determining cross-over
    cubic.curve(j,1) = cubic.b0(j) + (cubic.b1(j)*cubic.root_b(j)) + (cubic.b2(j)*(cubic.root_b(j)^2)) + (cubic.b3(j)*(cubic.root_b(j)^3));
end % end of cubic loop

% quartic equation
% a4*x^4 + a3*x^3 + a2*x^2 + a1x + a0 = 0
j = 1;
quartic.a3_tol = 10^-3;
for j = 1:param.npoints
    quartic.a0(j) = (tor.w/tran.w)^2;
    quartic.a1(j) = 0;
    quartic.a2(j) = -((tor.w^2)/(tran.w^2)) - (4*tran.xi*tor.xi*tor.w/tran.w) - 1 - (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).graph(j);
    quartic.a3(j) = 2*tor.xi*(tor.w/tran.w)*(param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).graph(j) + 2*tran.xi*(param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).graph(j);
    if quartic.a3(j) < quartic.a3_tol
        quartic.a3(j) = 0;
    end % end of tolerance loop
%    a4(j) = 1 + (rho*(B^4)/I1)*C22*A3st(j) + (rho*(B^2)/M1)*(rho*(B^4)/I1)*((C12^2)*A1st(j)*H2st(j) - C11*C22*A2st(j)*H1st(j))
    quartic.a4(j) = 1 + ((param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).graph(j)) + ((param.rho*(param.B^2)/param.M1)*(param.rho*(param.B^4)/param.I1)*((param.C12*param.C12*fder.A(1).graph(j)*fder.H(2).graph(j))-(param.C11*param.C22*fder.A(2).graph(j)*fder.H(1).graph(j))));
    % find the roots of the polynomial
    quartic.roots_a(:,1) = roots([quartic.a4(j) quartic.a3(j) quartic.a2(j) quartic.a1(j) quartic.a0(j)]);
    k = 1;
    h = 1;
    for k=1:4
        if quartic.roots_a(k,1)>0
            quartic.root_a_sq(h,j) = quartic.roots_a(k,1);
%             if k == 2
%                 root_a(h,j) = sqrt(root_a_sq(h,j))
%                 else
            quartic.root_a(h,j) = quartic.root_a_sq(h,j);
%            end
            h = h+1;
        end % end of positive root if statement
    end % end of positive-root-finding loop
    q = h-1;
    % assemble points for determining cross-over
    p = 1;
    for p = 1:q;
        quartic.curve(j,p) = quartic.a0(j) + (quartic.a1(j)*quartic.root_a(p,j)) + (quartic.a2(j)*(quartic.root_a(p,j)^2)) + (quartic.a3(j)*(quartic.root_a(p,j)^3)) + (quartic.a4(j)*(quartic.root_a(p,j)^4));
    end
    
end % end of quartic loop

% assemble points for determining cross-over
j = 1;
for j = 1:param.npoints  
    fder.K(1,j) = 2*pi/fder.U_NB.graph(j);
%        fder.K(j,1) = 2*pi/fder.U_NB.graph(j);
end

% plot roots, X, versus reduced frequency, K
figure
%plot(K,cubic,'-.r',K,quartic(:,1),'-g',K,quartic(:,2),':b','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency'),legend('Cubic','Quartic 1','Quartic 2')
plot(fder.K,cubic.root_b','-.r',fder.K,quartic.root_a(1,:)','-g',fder.K,quartic.root_a(2,:)','-','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency: X versus K'),legend('Cubic','Quartic 1','Quartic 2')
axis([0.5 pi 0.95 2.05])
grid on

% detailed plot of roots, X, versus reduced frequency, K
figure
%plot(K,cubic,'-.r',K,quartic(:,1),'-g',K,quartic(:,2),':b','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency'),legend('Cubic','Quartic 1','Quartic 2')
plot(fder.K,cubic.root_b','-.r',fder.K,quartic.root_a(1,:)','-g',fder.K,quartic.root_a(2,:)','-','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency: X versus K, Detailed Section'),legend('Cubic','Quartic 1')
%plot(K((npoints-1):npoints,:),root_b(1,(npoints-1):npoints)','-.r',K((npoints-1):npoints,:),root_a(1,(npoints-1):npoints)','-g','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency: X versus K, Detailed Section'),legend('Cubic','Quartic 1')
axis([0.52 0.64 1.6 1.68])
grid on

% Critical Flutter Speed, from numerical example (reworked)
% determined automatically !!
[fder.Kc.graph,fder.Xc.graph] = polyxpoly(fder.K,cubic.root_b',fder.K,quartic.root_a(1,:)'); 
fder.Uc.graph = param.B*tran.w*fder.Xc.graph/fder.Kc.graph;
fder.Uc.graph % 183.9912
fder.Uc_mph.graph = fder.Uc.graph*3600/5280;
fder.Uc_mph.graph %  125.4485

% determined manually from plot
% Xc = 1.625
% Kc = 0.555

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution of the Flutter Equations using theoretical flutter derivatives %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trial 1:  Using only the 6 points from the table %%%%%%%%%%%%%%%%%%%%%%%%
fder.K.graph = fder.K;
fder.K.TCF = fder.K.graph;
clear fder.K;

for j = 1:param.npoints
   % K = K_s(j); % value of reduced frequency to be used in iteration j
    z = fder.K.graph(j)/2; % here z is a function of k = B*w/(2U) not K!
    %% THEODORSEN CIRCULATORY FUNCTION C(k) %%
    % real part of Theodorsen Circulatory Function (TCF)
    TCF.F(j) = (besselj(1,z)*(besselj(1,z)+bessely(0,z))+bessely(1,z)*(bessely(1,z)-besselj(0,z)))/((besselj(1,z)+bessely(0,z))^2+(bessely(1,z)-besselj(0,z))^2);
    % imaginary part of Theodorsen Circulatory Function (TCF)
    TCF.G(j) = -(bessely(1,z)*bessely(0,z)+besselj(1,z)*besselj(0,z))/((besselj(1,z)+bessely(0,z))^2+(bessely(1,z)-besselj(0,z))^2);
    % Complex Theodorsen Circulatory Function
    TCF.C(j) = TCF.F(j) + i*TCF.G(j);

    % DEFINE FLUTTER DERIVATIVES/AERODYNAMIC COEFFICIENTS
    fder.H(1).TCF(j) = -pi*TCF.F(j)/fder.K.graph(j);
    fder.H(2).TCF(j) = -(pi/(4*fder.K.graph(j)))*(1+(4*TCF.G(j)/fder.K.graph(j))+TCF.F(j));
    fder.H(3).TCF(j) = -(pi/(2*(fder.K.graph(j)^2)))*((2*TCF.F(j))-(TCF.G(j)*fder.K.graph(j)/2));
    fder.H(4).TCF(j) = (pi/4)*(1+(4*TCF.G(j)/fder.K.graph(j)));
    
    fder.A(1).TCF(j) = pi*TCF.F(j)/(4*fder.K.graph(j));
    fder.A(2).TCF(j) = -(pi/(4*fder.K.graph(j)*fder.K.graph(j)))*((fder.K.graph(j)/4)-TCF.G(j)-(fder.K.graph(j)*TCF.F(j)/4));
    fder.A(3).TCF(j) = (pi/(4*fder.K.graph(j)*fder.K.graph(j)))*((fder.K.graph(j)*fder.K.graph(j)/32)+TCF.F(j)-(TCF.G(j)*fder.K.graph(j)/4));
    fder.A(4).TCF(j) = -(pi/(4*fder.K.graph(j)))*TCF.G(j);

end

% add (0,0) point for plotting
fder.A0(1).TCF = [0 fder.A(1).TCF];
fder.A0(2).TCF = [0 fder.A(2).TCF];
fder.A0(3).TCF = [0 fder.A(3).TCF];
fder.A0(4).TCF = [0 fder.A(4).TCF];
fder.H0(1).TCF = [0 fder.H(1).TCF];
fder.H0(2).TCF = [0 fder.H(2).TCF];
fder.H0(3).TCF = [0 fder.H(3).TCF];
fder.H0(4).TCF = [0 fder.H(4).TCF];

% Compare Tabulated and Theoretical Values of Aerodynamic Coefficients
pos=0;
fder.U_NB0.graph
figure
subplot(3,2,1);plot(fder.U_NB0.graph,fder.A0(1).graph,fder.U_NB0.graph,fder.A0(1).TCF,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('A_1*'),legend('Scanlan (1981)','Theoretical');axis([0 12 0 3]);grid on
subplot(3,2,2);plot(fder.U_NB0.graph,fder.H0(1).graph,fder.U_NB0.graph,fder.H0(1).TCF,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('-H_1*');axis([0 12 0 6]);grid on
subplot(3,2,3);plot(fder.U_NB0.graph,fder.A0(2).graph,fder.U_NB0.graph,fder.A0(2).TCF,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('A_2*');axis([0 12 -0.2 0.4]);grid on
subplot(3,2,4);plot(fder.U_NB0.graph,fder.H0(2).graph,fder.U_NB0.graph,fder.H0(2).TCF,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('H_2*');axis([0 12 -2 8]);grid on
subplot(3,2,5);plot(fder.U_NB0.graph,fder.A0(3).graph,fder.U_NB0.graph,fder.A0(3).TCF,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('A_3*');axis([0 12 0 3]);grid on
subplot(3,2,6);plot(fder.U_NB0.graph,fder.H0(3).graph,fder.U_NB0.graph,fder.H0(3).TCF,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('-H_3*');axis([0 12 0 8]);grid on

%% FIND ROOTS OF COMPLEX POLYNOMIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cubic equation
% b3*x^3 + b2*x^2 + b1x + b0 = 0
j = 1;
for j = 1:param.npoints
    cubic.b0t(j) = 2*tran.xi*((tor.w/tran.w)^2) + 2*tor.xi*tor.w/tran.w;
    cubic.b1t(j) = -(param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j)*((tor.w/tran.w)^2) - (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j);
    cubic.b2t(j) = -2*tor.xi*tor.w/tran.w - 2*tran.xi - 2*tran.xi*(param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j);
    cubic.b3t(j) = (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j) + (param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j) + (param.rho*(param.B^2)/param.M1)*(param.rho*(param.B^4)/param.I1)*(param.C11*param.C22*fder.H(1).TCF(j)*fder.A(3).TCF(j) - (param.C12^2)*fder.A(1).TCF(j)*fder.H(3).TCF(j));
    cubic.roots_bt(:,1) = roots([cubic.b3t(j) cubic.b2t(j) cubic.b1t(j) cubic.b0t(j)]);
    k = 1;
    for k=1:3
        if cubic.roots_bt(k,1)>0
            cubic.root_bt(j) = cubic.roots_bt(k,1);
        end % end of positive root if statement
    end % end of positive-root-finding loop
    % assemble points for determining cross-over
    cubic.curvet(j,1) = cubic.b0t(j) + (cubic.b1t(j)*cubic.root_bt(j)) + (cubic.b2t(j)*(cubic.root_bt(j)^2)) + (cubic.b3t(j)*(cubic.root_bt(j)^3));
end % end of cubic loop

% quartic equation
% a4*x^4 + a3*x^3 + a2*x^2 + a1x + a0 = 0
j = 1;
quartic.a3_tol = 10^-3;
for j = 1:param.npoints
    quartic.a0t(j) = (tor.w/tran.w)^2;
    quartic.a1t(j) = 0;
    quartic.a2t(j) = -((tor.w^2)/(tran.w^2)) - (4*tran.xi*tor.xi*tor.w/tran.w) - 1 - (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j);
    quartic.a3t(j) = 2*tor.xi*(tor.w/tran.w)*(param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j) + 2*tran.xi*(param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j);
    if quartic.a3t(j) < quartic.a3_tol
        quartic.a3t(j) = 0;
    end % end of tolerance loop
%    a4(j) = 1 + (rho*(B^4)/I1)*C22*A3st(j) + (rho*(B^2)/M1)*(rho*(B^4)/I1)*((C12^2)*A1st(j)*H2st(j) - C11*C22*A2st(j)*H1st(j))
    quartic.a4t(j) = 1 + ((param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j)) + ((param.rho*(param.B^2)/param.M1)*(param.rho*(param.B^4)/param.I1)*((param.C12*param.C12*fder.A(1).TCF(j)*fder.H(2).TCF(j))-(param.C11*param.C22*fder.A(2).TCF(j)*fder.H(1).TCF(j))));
    % find the roots of the polynomial
    quartic.roots_at(:,1) = roots([quartic.a4t(j) quartic.a3t(j) quartic.a2t(j) quartic.a1t(j) quartic.a0t(j)]);
    k = 1;
    h = 1;
    for k=1:4
        if quartic.roots_at(k,1)>0
            quartic.root_a_sqt(h,j) = quartic.roots_at(k,1);
%             if k == 2
%                 root_a(h,j) = sqrt(root_a_sq(h,j))
%                 else
            quartic.root_at(h,j) = quartic.root_a_sqt(h,j);
%            end
            h = h+1;
        end % end of positive root if statement
    end % end of positive-root-finding loop
    q = h-1;
    % assemble points for determining cross-over
    p = 1;
    for p = 1:q;
        quartic.curvet(j,p) = quartic.a0t(j) + (quartic.a1t(j)*quartic.root_at(p,j)) + (quartic.a2t(j)*(quartic.root_at(p,j)^2)) + (quartic.a3t(j)*(quartic.root_at(p,j)^3)) + (quartic.a4t(j)*(quartic.root_at(p,j)^4));
    end
end % end of quartic loop

% plot roots, X, versus reduced frequency, K
figure
%plot(K,cubic,'-.r',K,quartic(:,1),'-g',K,quartic(:,2),':b','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency'),legend('Cubic','Quartic 1','Quartic 2')
plot(fder.K.TCF,cubic.root_bt','-.r',fder.K.TCF,quartic.root_at(1,:)','-g',fder.K.TCF,quartic.root_at(2,:)','-','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency: X versus K'),legend('Cubic','Quartic 1','Quartic 2')
axis([0.5 pi 0.95 2.05])
grid on

% % detailed plot of roots, X, versus reduced frequency, K
% figure
% %plot(K,cubic,'-.r',K,quartic(:,1),'-g',K,quartic(:,2),':b','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency'),legend('Cubic','Quartic 1','Quartic 2')
% plot(K,root_bt','-.r',K,root_at(1,:)','-g',K,root_at(2,:)','-','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency: X versus K, Detailed Section'),legend('Cubic','Quartic 1')
% %plot(K((npoints-1):npoints,:),root_b(1,(npoints-1):npoints)','-.r',K((npoints-1):npoints,:),root_a(1,(npoints-1):npoints)','-g','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency: X versus K, Detailed Section'),legend('Cubic','Quartic 1')
% axis([0.52 0.64 1.6 1.68])
% grid on

% % Trial 2:  Evaluation for many points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify number of points for U_NB
param.points = 70;
fder.U_NB.TCF = linspace(1,param.points,param.points);
fder.K.TCF = 2*pi./fder.U_NB.TCF;

for j = 1:param.points
   % K = K_s(j); % value of reduced frequency to be used in iteration j
    z = fder.K.TCF(j)/2; % here z is a function of k = B*w/(2U) not K!
    %% THEODORSEN CIRCULATORY FUNCTION C(k) %%
    % real part of Theodorsen Circulatory Function (TCF)
    TCF.F(j) = (besselj(1,z)*(besselj(1,z)+bessely(0,z))+bessely(1,z)*(bessely(1,z)-besselj(0,z)))/((besselj(1,z)+bessely(0,z))^2+(bessely(1,z)-besselj(0,z))^2);
    % imaginary part of Theodorsen Circulatory Function (TCF)
    TCF.G(j) = -(bessely(1,z)*bessely(0,z)+besselj(1,z)*besselj(0,z))/((besselj(1,z)+bessely(0,z))^2+(bessely(1,z)-besselj(0,z))^2);
    % Complex Theodorsen Circulatory Function
    TCF.C(j) = TCF.F(j) + i*TCF.G(j);

    % DEFINE FLUTTER DERIVATIVES/AERODYNAMIC COEFFICIENTS
    fder.H(1).TCF(j) = -pi*TCF.F(j)/fder.K.TCF(j);
    fder.H(2).TCF(j) = -(pi/(4*fder.K.TCF(j)))*(1+(4*TCF.G(j)/fder.K.TCF(j))+TCF.F(j));
    fder.H(3).TCF(j) = -(pi/(2*(fder.K.TCF(j)^2)))*((2*TCF.F(j))-(TCF.G(j)*fder.K.TCF(j)/2));
    fder.H(4).TCF(j) = (pi/4)*(1+(4*TCF.G(j)/fder.K.TCF(j)));
    
    fder.A(1).TCF(j) = pi*TCF.F(j)/(4*fder.K.TCF(j));
    fder.A(2).TCF(j) = -(pi/(4*fder.K.TCF(j)*fder.K.TCF(j)))*((fder.K.TCF(j)/4)-TCF.G(j)-(fder.K.TCF(j)*TCF.F(j)/4));
    fder.A(3).TCF(j) = (pi/(4*fder.K.TCF(j)*fder.K.TCF(j)))*((fder.K.TCF(j)*fder.K.TCF(j)/32)+TCF.F(j)-(TCF.G(j)*fder.K.TCF(j)/4));
    fder.A(4).TCF(j) = -(pi/(4*fder.K.TCF(j)))*TCF.G(j);

end

%% FIND ROOTS OF COMPLEX POLYNOMIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cubic equation
% b3*x^3 + b2*x^2 + b1x + b0 = 0
j = 1;
for j = 1:param.points
    cubic.b0t(j) = 2*tran.xi*((tor.w/tran.w)^2) + 2*tor.xi*tor.w/tran.w;
    cubic.b1t(j) = -(param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j)*((tor.w/tran.w)^2) - (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j);
    cubic.b2t(j) = -2*tor.xi*tor.w/tran.w - 2*tran.xi - 2*tran.xi*(param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j);
    cubic.b3t(j) = (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j) + (param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j) + (param.rho*(param.B^2)/param.M1)*(param.rho*(param.B^4)/param.I1)*(param.C11*param.C22*fder.H(1).TCF(j)*fder.A(3).TCF(j) - (param.C12^2)*fder.A(1).TCF(j)*fder.H(3).TCF(j));
    cubic.roots_bt(:,1) = roots([cubic.b3t(j) cubic.b2t(j) cubic.b1t(j) cubic.b0t(j)]);
    k = 1;
    for k=1:3
        if cubic.roots_bt(k,1)>0
            cubic.root_bt(j) = cubic.roots_bt(k,1);
        end % end of positive root if statement
    end % end of positive-root-finding loop
    % assemble points for determining cross-over
    cubic.curvet(j,1) = cubic.b0t(j) + (cubic.b1t(j)*cubic.root_bt(j)) + (cubic.b2t(j)*(cubic.root_bt(j)^2)) + (cubic.b3t(j)*(cubic.root_bt(j)^3));
end % end of cubic loop

% quartic equation
% a4*x^4 + a3*x^3 + a2*x^2 + a1x + a0 = 0
j = 1;
quartic.a3_tol = 10^-3;
for j = 1:param.points
    quartic.a0t(j) = (tor.w/tran.w)^2;
    quartic.a1t(j) = 0;
    quartic.a2t(j) = -((tor.w^2)/(tran.w^2)) - (4*tran.xi*tor.xi*tor.w/tran.w) - 1 - (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j);
    quartic.a3t(j) = 2*tor.xi*(tor.w/tran.w)*(param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j) + 2*tran.xi*(param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j);
    if quartic.a3t(j) < quartic.a3_tol
        quartic.a3t(j) = 0;
    end % end of tolerance loop
%    a4(j) = 1 + (rho*(B^4)/I1)*C22*A3st(j) + (rho*(B^2)/M1)*(rho*(B^4)/I1)*((C12^2)*A1st(j)*H2st(j) - C11*C22*A2st(j)*H1st(j))
    quartic.a4t(j) = 1 + ((param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j)) + ((param.rho*(param.B^2)/param.M1)*(param.rho*(param.B^4)/param.I1)*((param.C12*param.C12*fder.A(1).TCF(j)*fder.H(2).TCF(j))-(param.C11*param.C22*fder.A(2).TCF(j)*fder.H(1).TCF(j))));
    % find the roots of the polynomial
    quartic.roots_at(:,1) = roots([quartic.a4t(j) quartic.a3t(j) quartic.a2t(j) quartic.a1t(j) quartic.a0t(j)]);
    k = 1;
    h = 1;
    for k=1:4
        if quartic.roots_at(k,1)>0
            quartic.root_a_sqt(h,j) = quartic.roots_at(k,1);
%             if k == 2
%                 root_a(h,j) = sqrt(root_a_sq(h,j))
%                 else
            quartic.root_at(h,j) = quartic.root_a_sqt(h,j);
%            end
            h = h+1;
        end % end of positive root if statement
    end % end of positive-root-finding loop
    q = h-1;
    % assemble points for determining cross-over
    p = 1;
    for p = 1:q;
        quartic.curvet(j,p) = quartic.a0t(j) + (quartic.a1t(j)*quartic.root_at(p,j)) + (quartic.a2t(j)*(quartic.root_at(p,j)^2)) + (quartic.a3t(j)*(quartic.root_at(p,j)^3)) + (quartic.a4t(j)*(quartic.root_at(p,j)^4));
    end
end % end of quartic loop

% plot roots, X, versus reduced frequency, K
figure
%plot(K,cubic,'-.r',K,quartic(:,1),'-g',K,quartic(:,2),':b','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency'),legend('Cubic','Quartic 1','Quartic 2')
plot(fder.K.TCF,cubic.root_bt','-.r',fder.K.TCF,quartic.root_at(1,:)','-g',fder.K.TCF,quartic.root_at(2,:)','-','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency: X versus K'),legend('Cubic','Quartic 1','Quartic 2')
axis([0 pi 0.95 2.05])
grid on

% detailed plot of roots, X, versus reduced frequency, K
figure
plot(fder.K.TCF,cubic.root_bt','-.r',fder.K.TCF,quartic.root_at(1,:)','-g','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency: X versus K, Detailed Section'),legend('Cubic','Quartic 1')
axis([0.3 0.50 1.4 1.5])
grid on

% Critical Flutter Speed, using 3 theorectical f.utter derivatices
[fder.Kc.TCF,fder.Xc.TCF] = polyxpoly(fder.K.TCF,cubic.root_bt',fder.K.TCF,quartic.root_at(1,:)'); 
fder.Uc.TCF = param.B*tran.w*fder.Xc.TCF/fder.Kc.TCF;
fder.Uc.TCF % 224.6868
fder.Uc_mph.TCF = fder.Uc.TCF*3600/5280;
fder.Uc_mph.TCF %  153.1955

% Comparison of Experimental Value with Theoretical Value %%%%%%%%%
change = 100*(fder.Uc.TCF-fder.Uc.graph)/fder.Uc.graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution Using ALL Flutter Derivatives (Theoretical Only) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot additional flutter derivatives
fder.U_NB0.TCF = [0 fder.U_NB.TCF];
fder.A0(4).TCF = [0 fder.A(4).TCF];
fder.H0(4).TCF = [0 fder.H(4).TCF];

figure
plot(fder.U_NB0.TCF,fder.A0(4).TCF,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('A_4*');axis([0 12 0 0.3]);grid on
figure
plot(fder.U_NB0.TCF,fder.H0(4).TCF,'LineWidth',2),xlabel('U/NB = 2*pi/K'),ylabel('H_4*');axis([0 12 -0.5 1]);grid on

% cubic equation
% b3*x^3 + b2*x^2 + b1x + b0 = 0
j = 1;
for j = 1:param.points
    cubic.b0t(j) = 2*tran.xi*((tor.w/tran.w)^2) + 2*tor.xi*tor.w/tran.w;
    cubic.b1t(j) = -(param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j)*((tor.w/tran.w)^2) - (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j);
    cubic.b2t4(j) = -2*tor.xi*tor.w/tran.w - 2*tran.xi - 2*tran.xi*(param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j)-2*tor.xi*(param.rho*(param.B^2)/param.M1)*(tor.w/tran.w)*param.C11*fder.H(4).TCF(j);
    cubic.b3t4(j) = (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j) + (param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j) + (param.rho*(param.B^2)/param.M1)*(param.rho*(param.B^4)/param.I1)*(param.C11*param.C22*(fder.H(1).TCF(j)*fder.A(3).TCF(j) + fder.A(2).TCF(j)*fder.H(4).TCF(j)) - (param.C12^2)*(fder.A(1).TCF(j)*fder.H(3).TCF(j) + fder.A(4).TCF(j)*fder.H(2).TCF(j)));
    cubic.roots_bt4(:,1) = roots([cubic.b3t4(j) cubic.b2t4(j) cubic.b1t(j) cubic.b0t(j)]);
    k = 1;
    for k=1:3
        if cubic.roots_bt4(k,1)>0
            cubic.root_bt4(j) = cubic.roots_bt4(k,1);
        end % end of positive root if statement
    end % end of positive-root-finding loop
    % assemble points for determining cross-over
    cubic.curvet4(j,1) = cubic.b0t(j) + (cubic.b1t(j)*cubic.root_bt4(j)) + (cubic.b2t4(j)*(cubic.root_bt4(j)^2)) + (cubic.b3t4(j)*(cubic.root_bt4(j)^3));
end % end of cubic loop

% quartic equation
% a4*x^4 + a3*x^3 + a2*x^2 + a1x + a0 = 0
j = 1;
quartic.a3_tol = 10^-3;
for j = 1:param.points
    quartic.a0t(j) = (tor.w/tran.w)^2;
    quartic.a1t(j) = 0;
    quartic.a2t4(j) = -((tor.w^2)/(tran.w^2)) - (4*tran.xi*tor.xi*tor.w/tran.w) - 1 - (param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j) - (param.rho*(param.B^2)/param.M1)*((tor.w/tran.w)^2)*param.C11*fder.H(4).TCF(j);
    quartic.a3t(j) = 2*tor.xi*(tor.w/tran.w)*(param.rho*(param.B^2)/param.M1)*param.C11*fder.H(1).TCF(j) + 2*tran.xi*(param.rho*(param.B^4)/param.I1)*param.C22*fder.A(2).TCF(j);
    if quartic.a3t(j) < quartic.a3_tol
        quartic.a3t(j) = 0;
    end % end of tolerance loop
    quartic.a4t4(j) = 1 + ((param.rho*(param.B^4)/param.I1)*param.C22*fder.A(3).TCF(j)) + ((param.rho*(param.B^2)/param.M1)*(param.rho*(param.B^4)/param.I1)*((param.C12*param.C12*(fder.A(1).TCF(j)*fder.H(2).TCF(j) - fder.A(4).TCF(j)*fder.H(3).TCF(j)))-(param.C11*param.C22*(fder.A(2).TCF(j)*fder.H(1).TCF(j)-fder.A(3).TCF(j)*fder.H(4).TCF(j)))));
    % find the roots of the polynomial
    quartic.roots_at4(:,1) = roots([quartic.a4t4(j) quartic.a3t(j) quartic.a2t4(j) quartic.a1t(j) quartic.a0t(j)]);
    k = 1;
    h = 1;
    for k=1:4
        if quartic.roots_at4(k,1)>0
            quartic.root_a_sqt4(h,j) = quartic.roots_at4(k,1);
%             if k == 2
%                 root_a(h,j) = sqrt(root_a_sq(h,j))
%                 else
            quartic.root_at4(h,j) = quartic.root_a_sqt4(h,j);
%            end
            h = h+1;
        end % end of positive root if statement
    end % end of positive-root-finding loop
    q = h-1;
    % assemble points for determining cross-over
    p = 1;
    for p = 1:q;
        quartic.curvet4(j,p) = quartic.a0t(j) + (quartic.a1t(j)*quartic.root_at4(p,j)) + (quartic.a2t4(j)*(quartic.root_at4(p,j)^2)) + (quartic.a3t(j)*(quartic.root_at4(p,j)^3)) + (quartic.a4t4(j)*(quartic.root_at4(p,j)^4));
    end
end % end of quartic loop

% plot roots, X, versus reduced frequency, K
figure
%plot(K,cubic,'-.r',K,quartic(:,1),'-g',K,quartic(:,2),':b','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency'),legend('Cubic','Quartic 1','Quartic 2')
plot(fder.K.TCF,cubic.root_bt4','-.r',fder.K.TCF,quartic.root_at4(1,:)','-g',fder.K.TCF,quartic.root_at4(2,:)','-','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency including A_4* and H_4*: X versus K'),legend('Cubic','Quartic 1','Quartic 2')
axis([0 pi 0.95 2.05])
grid on

% detailed plot of roots, X, versus reduced frequency, K
figure
plot(fder.K.TCF,cubic.root_bt4','-.r',fder.K.TCF,quartic.root_at4(1,:)','-g',fder.K.TCF,quartic.root_at4(2,:)','-','LineWidth',2),xlabel('Reduced Frequency, K = Bw/U'),ylabel('X, from cubic or quartic equation'),title('Determination of Flutter Frequency including A_4* and H_4*: X versus K, Detailed Section'),legend('Cubic','Quartic 1',pos)
axis([0.3 0.50 1.4 1.5])
grid on

% Critical Flutter Speed, from numerical example (reworked)
[fder.Kc.TCF4,fder.Xc.TCF4] = polyxpoly(fder.K.TCF,cubic.root_bt4',fder.K.TCF,quartic.root_at4(1,:)') 
fder.Uc.TCF4 = param.B*tran.w*fder.Xc.TCF4/fder.Kc.TCF4;
fder.Uc.TCF4 % 208.6644
fder.Uc_mph.TCF4 = fder.Uc.TCF4*3600/5280;
fder.Uc_mph.TCF4 % 142.2712

% Comparison of Experimental Value with Theoretical Value %%%%%%%%%
change4 = 100*(fder.Uc.TCF4-fder.Uc.graph)/fder.Uc.graph