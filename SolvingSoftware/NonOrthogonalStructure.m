%% stiffness matrix method
clc;
clear;
n = 5; % number of members
I=[1 1 1 1 1];
L = [5 5 3 5 2.5]; % length in m
A = [6 12 12 6 6]/100; % Area in m2
theta= [90 0 0 -90 -90]; % angle in degrees   
uu = 9; % Number of unrestrained degrees of freedom
ur = 9; % Number of restrained degrees of freedom
l1 = [10 1 11 2 12 3]; % Global labels for member 1
l2 = [1 4 2 5 3 6]; % Global labels for member 2
l3 = [4 7 5 8 6 9]; % Global labels for member 3
l4 = [4 13 5 14 6 15]; % Global labels for member 4
l5 = [7 16 8 17 9 18]; % Global labels for member 5
l6 = [11 14 10 13 12 15]; % Global labels for member 6
l7 = [17 14 23 13 24 15]; % Global labels for member 7
l8 = [1 4 2 5 3 6]; % Global labels for member 8
l9 = [1 4 2 5 3 6]; % Global labels for member 9
l10 = [1 4 2 5 3 6]; % Global labels for member 10
l11 = [1 4 2 5 3 6]; % Global labels for member 10
l12 = [1 4 2 5 3 6]; % Global labels for member 10
l13 = [1 4 2 5 3 6]; % Global labels for member 10
l14 = [1 4 2 5 3 6]; % Global labels for member 10
l15 = [1 4 2 5 3 6]; % Global labels for member 10
l16 = [1 4 2 5 3 6]; % Global labels for member 10
l17 = [1 4 2 5 3 6]; % Global labels for member 10
l18 = [1 4 2 5 3 6]; % Global labels for member 10
l19 = [1 4 2 5 3 6]; % Global labels for member 10
l20 = [1 4 2 5 3 6]; % Global labels for member 10
l= [l1; l2; l3; l4; l5; l6; l7; l8; l9; l10;l11;l12;l13;l14;l15;l16;l17;l18;l19;l20];
dof = uu + ur; % Degrees of freedom
Ktotal = zeros (dof);

%% Fem CALC
L1=L;
l1=l;
dof = uu+ur;
al=[0 0 0 0];
w(:,:,1) = [6 -25 5 0;al;al;al];
w(:,:,2) = [1 30 5 0;al;al;al];
w(:,:,3) = [2 20 3 1;al;al;al];
w(:,:,4) = [1 -20 5 0;al;al;al];
w(:,:,5) = [al;al;al;al];
w(:,:,6) = [al;al;al;al];
w(:,:,7) = [al;al;al;al];
w(:,:,8) = [al;al;al;al];
w(:,:,9) = [al;al;al;al];
w(:,:,10) = [al;al;al;al];
w(:,:,11) = [al;al;al;al];
w(:,:,12) = [al;al;al;al];
w(:,:,13) = [al;al;al;al];
w(:,:,14) = [al;al;al;al];
w(:,:,15) = [al;al;al;al];
w(:,:,16) = [al;al;al;al];
w(:,:,17) = [al;al;al;al];
w(:,:,18) = [al;al;al;al];
w(:,:,19) = [al;al;al;al];
w(:,:,20) = [al;al;al;al];
fem=zeros(n,6);
E=1E0;
I=I*E;
for i=1:n
for x=1:4
    r=w(x,1,i);
    if r==0
    break;
    elseif(r==1)
      F=w(x,2,i);L=w(x,3,i);
      fem(i,1)=fem(i,1)+F*L/8;
      fem(i,2)=fem(i,2)-F*L/8;
      fem(i,3)=fem(i,3)+F/2;
      fem(i,4)=fem(i,4)+F/2;
    elseif(r==2)
      F=w(x,2,i);a=w(x,3,i);b=w(x,4,i);l=a+b;
      fem(i,1)=fem(i,1)+F*a*b*b/(l*l);
      fem(i,2)=fem(i,2)-F*a*a*b/(l*l);
      fem(i,3)=fem(i,3)+F*b/l;
      fem(i,4)=fem(i,4)+F*a/l;
    elseif(r==3)
      W=w(x,2,i);L=w(x,3,i);
      fem(i,1)=fem(i,1)+W*L*L/12;
      fem(i,2)=fem(i,2)-W*L*L/12;
      fem(i,3)=fem(i,3)+W*L/2;
      fem(i,4)=fem(i,4)+W*L/2;
    elseif(r==4)
      W=w(x,2,i);L=w(x,3,i);
      fem(i,1)=fem(i,1)+5*W*L*L/96;
      fem(i,2)=fem(i,2)-5*W*L*L/96;
      fem(i,3)=fem(i,3)+W*L/4;
      fem(i,4)=fem(i,4)+W*L/4;
    elseif(r==5)
      W=w(x,2,i);L=w(x,3,i);
      fem(i,1)=fem(i,1)+W*L*L/20;
      fem(i,2)=fem(i,2)-W*L*L/30;
      fem(i,3)=fem(i,3)+W*L/3;
      fem(i,4)=fem(i,4)+W*L/6;
    elseif(r==6)
      W=w(x,2,i);L=w(x,3,i);
      fem(i,1)=fem(i,1)+W*L*L/30;
      fem(i,2)=fem(i,2)-W*L*L/20;
      fem(i,3)=fem(i,3)+W*L/6;
      fem(i,4)=fem(i,4)+W*L/3;
    end  
end
end
l= l1;L = L1;
disp("Label")
disp(l(1:n,:));
  
dof = uu+ur;
Ktotal = zeros (dof);
disp("FEM")
disp(fem);   



%% rotation coefficients for each member
rc1 = 4.*I./L;
rc2 = 2.*I./L;
rc3 = A./L;
cx = cosd(theta);
cy = sind(theta);


%% stiffness matrix 6 by 6
fembar = zeros(6,1,n);
Tt = zeros(6,6,n);
Kge = zeros(6,6,n);
for i = 1:n
        Knew = zeros (dof);
         k1 = [rc1(i); rc2(i); (rc1(i)+rc2(i))/L(i);
(-(rc1(i)+rc2(i))/L(i)); 0; 0];
       k2 = [rc2(i); rc1(i); (rc1(i)+rc2(i))/L(i);
(-(rc1(i)+rc2(i))/L(i)); 0; 0;];
k3 = [(rc1(i)+rc2(i))/L(i); (rc1(i)+rc2(i))/L(i);
(2*(rc1(i)+rc2(i))/(L(i)^2)); (-2*(rc1(i)+rc2(i))/(L(i)^2));
0; 0;];
        k4 = -k3;
        k5 = [0; 0; 0; 0; rc3(i); -rc3(i)];
        k6 = [0; 0; 0; 0; -rc3(i); rc3(i)];
        K = [k1 k2 k3 k4 k5 k6];
        fprintf ('Member Number =');
        disp (i);
        fprintf ('Local Stiffness matrix of member, [K] = \n');
        disp (K);
        T1 = [1; 0; 0; 0; 0; 0];
        T2 = [0; 1; 0; 0; 0; 0];
        T3 = [0; 0; cx(i); 0; cy(i); 0];
        T4 = [0; 0; 0; cx(i); 0; cy(i)];
        T5 = [0; 0; -cy(i); 0; cx(i); 0];
        T6 = [0; 0; 0; -cy(i); 0; cx(i)];
        T = [T1 T2 T3 T4 T5 T6];
        fprintf ('Tranformation matrix of member, [T] = \n');
        disp (T);
        Ttr = T';
        fprintf ('Tranformation matrix Transpose, [T] = \n');
        disp (Ttr);
        Kg = Ttr*K*T;
        fprintf ('Global Matrix, [K global] = \n');
        disp (Kg);
        for p = 1:6
            for q = 1:6
                Knew((l(i,p)),(l(i,q))) =Kg(p,q);
            end
        end
        Ktotal = Ktotal + Knew;
        Tt(:,:,i) = T;
        Kge(:,:,i) = Kg;
        
        fembar(:,1,:) = Tt(:,:,i)*fem';
        
end
fprintf ('Stiffness Matrix of complete structure, [Ktotal] = \n');
disp (Ktotal);
Kunr = zeros(uu);
for x=1:uu
        for y=1:uu
            Kunr(x,y)= Ktotal(x,y);
        end
end
fprintf ('Unrestrained Stiffness sub-matrix, [Kuu] = \n');
disp (Kunr);
KuuInv= inv(Kunr);
fprintf ('Inverse of Unrestrained Stiffness sub-matrix,[KuuInverse] = \n');
disp (KuuInv);
%% JLCALC
jl=[0; -10; 0; 0; -20; 0; 0; 0; 0 ;0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
c=1;
for i=1:n
for j=1:4
    jl(l(i,j))=jl(l(i,j))-fem(i,j);
end
end
L = L1;
jlu = jl(1:uu,1); % load vector in unrestrained dof
delu = KuuInv*jlu;
disp("JL");
disp(jl(1:dof,:));

disp("JLU");
disp(jlu);
fprintf ('Unrestrained displacements, [DelU] = \n');
disp (delu');
delr = zeros (ur,1);
del = zeros (dof,1);
del = [delu; delr];
deli= zeros (6,1);
delbar=zeros(6,1,n);
mbar = zeros(6,1,n);
for i = 1:n
    for p = 1:6
        deli(p,1) = del((l(i,p)),1) ;
    end
    delbar(:,:,i) = deli;
    mbar(:,:,i)= (Kge(:,:,i) * delbar(:,:,i)) + fembar (:,:,i);
    fprintf ('Member Number =');
    disp (i);
    fprintf('Global displacement matrix [DeltaBar] = \n')
    disp (delbar(:,:,i)');
    fprintf ('Global End moment matrix [MBar] = \n');
    disp (mbar(:,:,i)');  
       
end
%% check
for i = 1:dof
    jf(i,1)=0;
    for j = 1:n
        %disp(class(find(l(j,:)== i)));
        if isempty(find(l(j,:)==i))==0
            jf(i,1) = jf(i,1)+ mbar(find(l(j,:)==i),1,j);
        end
   
        
    end
end
disp('Joint Forces')
disp(jf')