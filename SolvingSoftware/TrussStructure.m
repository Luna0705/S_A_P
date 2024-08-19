%% stiffness matrix method for planar truss
% Input
clc;
clear;
n = 9; % number of members till 20
L = [3 3 3 5 5 5 5 sqrt(34) sqrt(34)]; % length in m
A = [1 1 1 0.8 0.8 0.8 0.8 1.5 1.5]; % Area in m2
theta= [90 90 90 0 0 0 0 45 -45 ]; % angle in degrees
uu = 8; % Number of unrestrained degrees-of-freedom
ur = 4; % Number of restrained degrees-of-freedom
l1 = [9 1 10 2]; % Global labels for member 1
l2 = [7 3 8 4]; % Global labels for member 2
l3 = [11 5 12 6]; % Global labels for member 3
l4 = [1 3 2 4]; % Global labels for member 4
l5 = [3 5 4 6];% Global labels for member 5
l6 = [9 7 10 8];% Global labels for member 6
l7 = [7 11 8 12];% Global labels for member 7
l8 = [9 3 10 4];% Global labels for member 8
l9 = [3 11 4 12];% Global labels for member 9

l10 = [11 15 12 16];% Global labels for member 10
l11 = [1 9 2 10];% Global labels for member 10
l12 = [9 5 10 6];% Global labels for member 10
l13 = [3 11  4 12];% Global labels for member 10
l14 = [11 7 12 8];% Global labels for member 10
l15 = [11 7 12 8];% Global labels for member 10
l16 = [11 7 12 8];% Global labels for member 10
l17 = [11 7 12 8];% Global labels for member 10
l18 = [11 7 12 8];% Global labels for member 10
l19 = [11 7 12 8];% Global labels for member 10
l20 = [11 7 12 8];% Global labels for member 10
l= [l1; l2; l3; l4; l5; l6; l7; l8; l9; l10;l11;l12;l13;l14;l15;l16;l17;l18;l19;l20];

dof = uu+ur; %Degree of Freedom

uul = zeros(1,uu);
for i=1:uu
    uul(1,i) = i;
end

url = zeros(1,ur);
for i=uu+1:ur+uu
    url(1,i) = i;
end

disp("Label");
disp(l(1:n,:));
Ktotal = zeros(dof);

fem=zeros(n,4);

%% Calc
Tt = zeros(4,4*n);

rc3 = A./L;
cx = cosd(theta);
cy = sind(theta);
kg = zeros(4, 4*n);

for i = 1:n
    Knew = zeros (dof);
    k1 = [0; 0; 0; 0];
    k2 = [0; 0; 0; 0];
    k3 = [0; 0; rc3(i); -rc3(i)];
    k4 = [0; 0; -rc3(i); rc3(i)];
    K = [k1 k2 k3 k4];
    fprintf( 'Member Number =');
    disp (i);
    fprintf( 'Local Stiffness Matrix of Member, [K] =\n');
    disp(K);
    T1 = [cx(i); 0; cy(i); 0];
    T2 = [0; cx(i); 0; cy(i)];
    T3 = [-cy(i); 0; cx(i); 0];
    T4 = [0; -cy(i); 0; cx(i)];
    T = [T1 T2 T3 T4];
    fprintf ('Transformation matrix of member, [T] =\n');
    disp (T);
    Ttr = T';
    fprintf('Transpose Transformtion matrix of member, [T] =\n');
    disp (Ttr);
    Kg = Ttr*K*T;
    fprintf( 'Global Matrix, [K global] = \n');
    disp (Kg);
    for p = 1:4
        for q = 1:4
            Knew ((l(i,p)),(l(i,q))) = Kg(p,q);
        end
    end
    Ktotal = Ktotal + Knew;
    
      Tt(1:4,4*i-3:4*i)= T;
      kg(1:4,4*i-3:4*i)= Kg;
      fembar(1:4,i)    = T'*fem(i,1:4)';
                    
end

fprintf ('Stiffness Matrix Of complete structure, [Ktotal] =\n');
disp (Ktotal);
Kunr = zeros(uu);
for x=1:uu
    for y=1:uu
        Kunr(x,y) = Ktotal(x,y);
    end
end
fprintf ('Unrestrained Stifness sub-matrix, [Kuu] =\n');
disp (Kunr);
KuuInv= inv(Kunr);
fprintf( 'Inverse of Unrestrained Stiffness Sub-Matrix, [KuuInverse] =\n');
disp (KuuInv);



%% Joint Load Vector
jl=[-20; 50; -80; 0; -20; 0; -50; 0; 0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]
jlu = jl(1:uu,1);%input the values of unrestrained part of joint load 
delu = KuuInv*jlu;
disp("JL");
disp(jl(1:dof,:));
disp("JLU");
disp(jlu);



%% CALC2
fprintf( 'Unrestrained displacements, [Delu] =\n');
disp (delu);
delr = zeros(ur,1);
del = zeros(dof,1);
del = [delu;delr];
deli = zeros(4 ,1);

for i = 1:n
    for p = 1:4
        deli(p,1) = del((l(i,p)),1);
    end
    delbar(1:4,i) = deli;
    mbar(1:4,i)   = kg(1:4,4*i-3:4*i)*delbar(1:4,i) + fembar(1:4,i);
    fprintf('Member Number = ');
    disp(i);
    fprintf('Global displacement of matrix [Deltabar] = \n');
    disp(delbar(1:4,i));
    fprintf('Global End moment matrix [MBar] = \n');
    disp(mbar(1:4,i));
end
dof=uu+ur;
for i = 1:dof
    jf(i,1)=0;
    for j = 1:n
        %disp(class(find(l(j,:)== i)));
        if isempty(find(l(j,:)==i))==0
            jf(i,1) = jf(i,1)+ mbar(find(l(j,:)==i),j);
        end
   
        
    end
end
fprintf("Joint forces = ")
disp(jf')                      