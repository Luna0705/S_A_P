 clc;
 clear;
 n = 3; % number of members
 EI = [1 1 1 ]; %Flexural rigidity
 EIy = EI;
 EIz = EI;
 GI = [0.25 0.25 0.25].*EI; %Torsional constant
 EA = [0.25 0.25 0.25].*EI; %Axial rigidity
 L = [3 3 3]; % length in m
 nj = n+1; % Number of Joints
 codm = [0 0 0; 3 0 0; 3 0 -3; 3 -3 -3]; %Coordinate wrt X,Y.Z: size=nj;
 dc = [1 0 0; 0 0 -1; 0 -1 0]; % Direction cosines for each member 
 tytr = [1 1 2]; % Type of transformation fo each member 
 psi = [0 0 90]; % Psi angle in degrees for each member
% C matrix
 c=zeros(3,3,n);
 c(:,:,1) = [1 0 0; 0 1 0; 0 0 1]; % C matrix for member 1
c(:,:,2) = [0 0 -1; 0 1 0; 1 0 0]; % C matrix for member 2
c(:,:,3) = [0 -1 0; -1 0 0; 0 0 -1]; % C matrix for member 3
% c(:,:,4) = [1 0 0; 0 1 0; 0 0 1]; % C matrix for member 4
% c(:,:,5) = [1 0 0; 0 1 0; 0 0 1]; % C matrix for member 5uu = 12; % Number of unrestrained Degrees-of-freedom
uu = 12; % Number of unrestrained Degrees-of-freedom
ur = 24; % Number of restrained Degrees-of-freedom uul = [1 2 3 4 5 6 7 8 9 10 11 12]; % global labels of unrestrained dof
 url = [13 14 15 16 17 18 19 20 21 22 23 24]; % global labels of restrained dof
 l1 = [13 14 15 16 17 18 1 2 3 4 5 6]; % Global labels formember 1
l2 = [1 2 3 4 5 6 19 20 21 22 23 24]; % Global labels formember 2
l3 = [25 26 27 28 29 30 7 8 9 10 11 12]; % Global labels formember 3
l4 = [7 8 9 10 11 12 31 32 33 34 35 36]; % Global labels formember 4
l5 = [7 8 9 10 11 12 1 2 3 4 5 6]; % Global labels for member5
l= [l1; l2; l3; l4; l5];
 dof = uu + ur; % Degrees-of-freedom
 Ktotal = zeros (dof);
 T = zeros(12,12,n);
 for i = 1:3
    for j = 1:3
        for k= 1:n

            T(i,j,k)=c(i,j,k);
            T(i+3,j+3,k)=c(i,j,k);
            T(i+6,j+6,k)=c(i,j,k);
            T(i+9,j+9,k)=c(i,j,k);

        end    
    end   
end
%% Fem Calc
l1= l;L1 = L;
al=[0 0 0 0 0];
w(:,:,1) = [1 40 4 0 1;al;al;al;al];
w(:,:,2) = [1 40 2 0 2;al;al;al;al];
w(:,:,3) = [al;al;al;al;al];
w(:,:,4) = [5 40 2 0 1;al;al;al;al];
w(:,:,5) = [al;al;al;al;al];
w(:,:,6) = [al;al;al;al;al];
w(:,:,7) = [al;al;al;al;al];
w(:,:,8) = [al;al;al;al;al];
w(:,:,9) = [al;al;al;al;al];
w(:,:,10) = [al;al;al;al;al];
w(:,:,11) = [al;al;al;al;al];
w(:,:,12) = [al;al;al;al;al];
w(:,:,13) = [al;al;al;al;al];
w(:,:,14) = [al;al;al;al;al];
w(:,:,15) = [al;al;al;al;al];
w(:,:,16) = [al;al;al;al;al];
w(:,:,17) = [al;al;al;al;al];
w(:,:,18) = [al;al;al;al;al];
w(:,:,19) = [al;al;al;al;al];
w(:,:,20) = [al;al;al;al;al];
%1-for y-axis , rest is for z-axis
fem=zeros(n,12);
for i=1:n
for x=1:5
    r=w(x,1,i);
    m=w(x,5,i);
    if r==0
    break;
    elseif(r==1)
      F=w(x,2,i);L=w(x,3,i);
      if m==1
      fem(i,6)=fem(i,6)+F*L/8;
      fem(i,12)=fem(i,12)-F*L/8;
      fem(i,2)=fem(i,2)+F/2;
      fem(i,8)=fem(i,8)+F/2;
      else
      fem(i,5)=fem(i,5)+F*L/8;
      fem(i,11)=fem(i,11)-F*L/8;
      fem(i,3)=fem(i,3)+F/2;
      fem(i,9)=fem(i,9)+F/2;
      end
    elseif(r==2)
        if m==1
      F=w(x,2,i);a=w(x,3,i);b=w(x,4,i);l=a+b;
      fem(i,6)=fem(i,1)+F*a*b*b/(l*l);
      fem(i,12)=fem(i,2)-F*a*a*b/(l*l);
      fem(i,2)=fem(i,3)+F*b/l;
      fem(i,8)=fem(i,4)+F*a/l;
        else
            F=w(x,2,i);a=w(x,3,i);b=w(x,4,i);l=a+b;
      fem(i,5)=fem(i,1)+F*a*b*b/(l*l);
      fem(i,11)=fem(i,2)-F*a*a*b/(l*l);
      fem(i,3)=fem(i,3)+F*b/l;
      fem(i,9)=fem(i,4)+F*a/l;
        end
    elseif(r==3)
        if m==1
      W=w(x,2,i);L=w(x,3,i);
      fem(i,6)=fem(i,1)+W*L*L/12;
      fem(i,12)=fem(i,2)-W*L*L/12;
      fem(i,2)=fem(i,3)+W*L/2;
      fem(i,8)=fem(i,4)+W*L/2;
        else
      W=w(x,2,i);L=w(x,3,i);
      fem(i,5)=fem(i,1)+W*L*L/12;
      fem(i,11)=fem(i,2)-W*L*L/12;
      fem(i,3)=fem(i,3)+W*L/2;
      fem(i,9)=fem(i,4)+W*L/2;
        end
    elseif(r==4)
        if m==1
      W=w(x,2,i);L=w(x,3,i);
      fem(i,6)=fem(i,1)+5*W*L*L/96;
      fem(i,12)=fem(i,2)-5*W*L*L/96;
      fem(i,2)=fem(i,3)+W*L/4;
      fem(i,8)=fem(i,4)+W*L/4;
        else
      W=w(x,2,i);L=w(x,3,i);
      fem(i,5)=fem(i,1)+5*W*L*L/96;
      fem(i,11)=fem(i,2)-5*W*L*L/96;
      fem(i,3)=fem(i,3)+W*L/4;
      fem(i,9)=fem(i,4)+W*L/4;
        end
    elseif(r==5)
        if m==1
     W=w(x,2,i);L=w(x,3,i);
      fem(i,6)=fem(i,1)+W*L*L/20;
      fem(i,12)=fem(i,2)-W*L*L/30;
      fem(i,2)=fem(i,3)+W*L/3;
      fem(i,8)=fem(i,4)+W*L/6;
        else
     W=w(x,2,i);L=w(x,3,i);
      fem(i,5)=fem(i,1)+W*L*L/20;
      fem(i,11)=fem(i,2)-W*L*L/30;
      fem(i,3)=fem(i,3)+W*L/3;
      fem(i,9)=fem(i,4)+W*L/6;
        end
      
    elseif(r==6)
        if m==1
     W=w(x,2,i);L=w(x,3,i);
      fem(i,6)=fem(i,1)+W*L*L/30;
      fem(i,12)=fem(i,2)-W*L*L/20;
      fem(i,2)=fem(i,3)+W*L/6;
      fem(i,8)=fem(i,4)+W*L/3;
        else
     W=w(x,2,i);L=w(x,3,i);
      fem(i,5)=fem(i,1)+W*L*L/30;
      fem(i,11)=fem(i,2)-W*L*L/20;
      fem(i,3)=fem(i,3)+W*L/6;
      fem(i,9)=fem(i,4)+W*L/3;
        end
      
    end  
end
end

disp('Fem')
disp(fem);

Ktotal = zeros (dof);
l= l1;L = L1;

%% Getting Type of transformation and Psi angle
 for i = 1:n
  if tytr(i) ==1
    fprintf ('Member Number =');
    disp (i);
    fprintf ('Type of transformation is Y-Z-X \n');
  else
     fprintf ('Member Number =');
    disp (i);
    fprintf ('Type of transformation is Z-Y-X \n');
  end
  fprintf ('Psi angle=');
  disp (psi(i));
 end
%% Stiffness coefficients for each member
 sc1 = EA./L;
 sc2 = 6*EIz./(L.^2);
 sc3 = 6*EIy./(L.^2);
 sc4 = GI./L;
 sc5 = 2*EIy./L;
 sc6 = 12*EIz./(L.^3);
 sc7 = 12*EIy./(L.^3);
 sc8 = 2*EIz./L;
fem1=fem(1,:);
%% stiffness matrix 6 by 6
Kg=zeros(12,12,n);
Tt=zeros(12,12,n);
 Fembar=zeros(12,n);
 for i = 1:n
    Knew = zeros (dof);
    k1 = [sc1(i); 0; 0; 0; 0; 0; -sc1(i); 0; 0; 0; 0; 0];
    k2 = [0; sc6(i); 0; 0; 0; sc2(i); 0; -sc6(i); 0; 0; 0; sc2(i)];
    k3 = [0; 0; sc7(i); 0; -sc3(i); 0; 0; 0; -sc7(i); 0; -sc3(i); 0];
    k4 = [0; 0; 0; sc4(i); 0; 0; 0; 0; 0; -sc4(i); 0; 0];
    k5 = [0; 0; -sc3(i); 0; (2*sc5(i)); 0; 0; 0; sc3(i); 0; sc5(i); 0];
    k6 = [0; sc2(i); 0; 0; 0; (2*sc8(i)); 0; -sc2(i); 0; 0; 0; sc8(i)];
    k7 = -k1;
    k8 = -k2;
    k9 = -k3;
    k10 = -k4;
    k11 = [0; 0; -sc3(i); 0; sc5(i); 0; 0; 0; sc3(i); 0; (2*sc5(i)); 0];
    k12 = [0; sc2(i); 0; 0; 0; sc8(i); 0; -sc2(i); 0; 0; 0; (2*sc8(i))];
    K = [k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12];
    fprintf ('Member Number =');
    disp (i);
    fprintf ('Local Stiffness matrix of member, [K] = \n');
    disp (K);
    TK = T(:,:,i);
    Ttr = TK';
    Kg(:,:,i) = Ttr*K*TK;
    
    fprintf ('Transformation matrix, [T] = \n');
    disp (T(:,:,i));
    fprintf ('Global Matrix, [K global] = \n');
    disp (Kg(:,:,i));
    Kpid=Kg(:,:,i);
    for p = 1:12
      for q = 1:12
        Knew((l(i,p)),(l(i,q))) =Kpid(p,q);
      end
    end
    Ktotal = Ktotal + Knew;
      Fembar(:,i)= T(:,:,i)'*fem(i,:)';
    
 end
 fprintf ('Stiffness Matrix of complete structure, [Ktotal] = \n');
 disp (Ktotal);
 Kunr = zeros(12);
 for x=1:uu
 for y=1:uu
    Kunr(x,y)= Ktotal(x,y);
 end
 end
 fprintf ('Unrestrained Stiffness sub-matrix, [Kuu] = \n');
 disp (Kunr);
 KuuInv= inv(Kunr);
 fprintf ('Inverse of Unrestrained Stiffness sub-matrix, [KuuInverse] = \n');
 disp (KuuInv);
 % Kg1=Kg(:,:,1);
 % Kg2=Kg(:,:,2);
 % Kg3=Kg(:,:,3);
 % Kg4=Kg(:,:,4);
 % Kg5=Kg(:,:,5);
 % T5=T(:,:,5);
%% Creation of joint load vector
 mbar=zeros(12,n);
 jl=[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
c=1;
for i=1:n
for j=1:12
    jl(l(i,j))=jl(l(i,j))-fem(i,j);
end
end
 jlu = jl(1:uu,1); % load vector in unrestrained dof 
 delu = KuuInv*jlu;
 fprintf ('Joint Load vector, [Jl] = \n');
 disp (jl);
 fprintf ('Unrestrained displacements, [DelU] = \n');
 disp (delu);
 delr = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
 del = zeros (dof,1);
 del = [delu; delr];
 deli= zeros (12,1);
 for i = 1:n
  for p = 1:12
    deli(p,1) = del((l(i,p)),1) ;
  end
   
      delbar = deli;
      mbar(:,i)= (Kg(:,:,i) * delbar)+Fembar(:,i);
      fprintf ('Member Number =');
      disp (i);
      fprintf ('Global displacement matrix [DeltaBar] = \n');
      disp (delbar);
      fprintf ('Global End moment matrix [MBar] = \n');
      disp (mbar(:,i)');
 end

%% check
 jf = zeros(dof,1);
 M=mbar';
 for a=1:n
   for b=1:12 % size of k matrix
    d = l(a,b);
    jfnew = zeros(dof,1);
    jfnew(d,1)=M(a,b);
    jf=jf+jfnew;
  end
 end
 fprintf ('Joint forces = \n');
 disp (jf);
