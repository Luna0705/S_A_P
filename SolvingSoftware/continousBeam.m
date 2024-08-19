%% stiffness matrix method
% Input
clc;
clear;
n = 3; % number of members
E=1E0;
uul=[];url=[];
I=[1 1 1];
L=[3 3 4];
I=I*E;
uu = 3; % Number of unrestrained degrees of freedom
ur = 5; % Number of restrained degrees of freedom
L1=L;
al=[0 0 0 0];
w(:,:,1) = [3 40 4 0;al;al;al];
w(:,:,2) = [3 40 2 0;al;al;al];
w(:,:,3) = [al;al;al;al];
w(:,:,4) = [al;al;al;al];
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
fem=zeros(n,4);

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
      F=w(x,2,i);a=w(x,3,i);b=w(x,4,i);L=a+b;
      fem(i,1)=fem(i,1)+F*a*b*b/(L*L);
      fem(i,2)=fem(i,2)-F*a*a*b/(L*L);
      fem(i,3)=fem(i,3)+F*b/L;
      fem(i,4)=fem(i,4)+F*a/L;
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
L=L1;
% fem(1,:)=[0 0 0 0]; % Local Fixed end moments of member 1

for i = 1:uu
    uul = [uul i];
end
url = [];
for i = 1:ur
    url = [url i+uu];
end

l1 = [1 5 4 7]; % Global labels for member 1
l2 = [1 2 6 10]; % Global labels for member 2
l3 = [2 8 4 9]; % Global labels for member 3
l4 = [2 3 10 12]; % Global labels for member 4
l5 = [3 11 4 13]; % Global labels for member 5
l6 = [2 8 4 9]; % Global labels for member 6
l7 = [2 3 10 12]; % Global labels for member 7
l8 = [3 11 4 13]; % Global labels for member 8
l9 = [3 11 4 13]; % Global labels for member 9
l10 = [3 11 4 13]; % Global labels for member 10
l11 = [3 11 4 13]; % Global labels for member 10
l12 = [3 11 4 13]; % Global labels for member 10
l13 = [3 11 4 13]; % Global labels for member 10
l14 = [3 11 4 13]; % Global labels for member 10
l15 = [3 11 4 13]; % Global labels for member 10
l16 = [3 11 4 13]; % Global labels for member 10
l17 = [3 11 4 13]; % Global labels for member 10
l18 = [3 11 4 13]; % Global labels for member 10
l19 = [3 11 4 13]; % Global labels for member 10
l20 = [3 11 4 13]; % Global labels for member 10
l= [l1; l2; l3; l4; l5; l6; l7; l8; l9; l10;l11;l12;l13;l14;l15;l16;l17;l18;l19;l20];
disp("Label")
disp(l(1:n,:));
  
dof = uu+ur;
Ktotal = zeros (dof);
disp("FEM")
disp(fem);   

%% Creation of joint load vector
% jl must be a column matrix
jl=[0; 0; 0; 0; 0; 0; 0; 0; 0 ;0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
for i = 1:dof
    jl(i,1)=0;
    for j = 1:n
        %disp(class(find(l(j,:)== i)));
        if isempty(find(l(j,:)==i))==0
            jl(i,1) = jl(i,1)- fem(j,find(l(j,:)==i));
        end
   
        
    end
end
disp("JL");
disp(jl(1:dof,:));
jlu = jl(1:uu,1); % load vector in unrestrained dof
disp("JLU");
disp(jlu);

%%rotation coeffecients for each member
rc1 = 4.*I./L;
rc2 = 2.*I./L;

%%stiffness matrix 4 by 4 (axial deformation neglected)
Kg = zeros(4,4,n);
for i = 1:n
    Knew = zeros (dof);
    k1 = [rc1(i); rc2(i); (rc1(i)+rc2(i))/L(i);
    (-(rc1(i)+rc2(i))/L(i))];
    k2 = [rc2(i); rc1(i); (rc1(i)+rc2(i))/L(i);
    (-(rc1(i)+rc2(i))/L(i))];
    k3 = [(rc1(i)+rc2(i))/L(i); (rc1(i)+rc2(i))/L(i);(2*(rc1(i)+rc2(i))/(L(i)^2)); (-2*(rc1(i)+rc2(i))/(L(i)^2))];
    k4 = -k3;
    K = [k1 k2 k3 k4];
    fprintf ('Member Number =');
    disp (i);
    fprintf ('Local Stiffness matrix of member, [K] = \n');
    disp (K);
    for p = 1:4
     for q = 1:4
       Knew((l(i,p)),(l(i,q))) =K(p,q);
     end
    end
    Ktotal = Ktotal + Knew;
    Kg(:,:,i)=K;
    
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
%% Calculation of displacements
delu = KuuInv*jlu;
fprintf ('Joint Load vector, [Jl] = \n');
disp (jl(1:dof,:)');
fprintf ('Unrestrained displacements, [DelU] = \n');
disp (delu');
delr = zeros (ur,1);
del = [delu; delr];
deli= zeros (4,1);


delbar=zeros(4,1,n);
mbar = zeros(4,1,n);
for i = 1:n
    for p = 1:4
        deli(p,1)= del((l(i,p)),1) ;
    end
    delbar(:,:,i) = deli;
    mbar(:,:,i)= (Kg(:,:,i) * delbar(:,:,i)) +fem(i,:)';
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
fprintf("Joint forces = ")
disp(jf')
        

    
