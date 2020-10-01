%%%%%%%%   lammps data file  %%%%%%%%%%
clc
clear
close all

%%%%%%%%   define parameters %%%%%%%%%%
mass=26010.98564;     %%原子量
r1=30;   %%%第1層r1
dr=1;   %%%層之間的厚度差
dn=2*dr; %%%層之間的顆數差
r2=r1-dr;  %%%第2層r2
l=2;    %%l為層數
nt=1;   %%%nt為atom種類
bt=1;   %%%bt為bond種類; 
nz1=60;   %%第一層particle數目
nz2=nz1-dn; %%第二層particle數目
nz=nz1+nz2; %%總層particle數目
hz=200;  %%hz為垂直高度原子數目
Ziv=3;     %%上下原子間隔差
h=Ziv*hz;    %%圓柱高度
ni=(nz1+nz2)*hz; %%atoms數目
l=2;      %%總層數
theta11=0;
theta12=2*pi/nz1;
theta13=theta11:theta12:2*pi-theta12;
theta21=0;
theta22=2*pi/nz2;
theta23=theta21:theta22:2*pi-theta22;
nb=l*nz*(hz-1)+ni; %%nb為總bond數目 
nangles=0;      %%angles數目
ndihedrals=0;   %%dihedrals數目
nimpropers=0;   %%impropers數目
anglet=0;     %%angles type
dihedralt=0;  %%dihedral type
impropert=0;  %%improper type



%%%%%%%%%%   box 大小 %%%%%%%%
xo=-1000;
xi=1000;
yo=-1000;
yi=1000;
ZO=-1000;
zi=1000;

%%%%%%%atom cooridnate%%%%%%%
%%%%%%%     X1,Y1建立   %%%%%%%
X1ii = [r1*cos(theta13)];
Y1ii = [r1*sin(theta13)];  
for j=1:1:nz1*hz
  g = rem(j,nz1); 
  if  g~=0
    X1i(j) = X1ii(1,g);
    Y1i(j) = Y1ii(1,g);
  else
    X1i(j) = X1ii(1,nz1);
    Y1i(j) = Y1ii(1,nz1);
  end
end
X1=X1i.';
Y1=Y1i.';
%%%%%%%     X2,Y2建立   %%%%%%%
X2ii = [r2*cos(theta23)];
Y2ii = [r2*sin(theta23)];  
for j=1:1:nz2*hz
  g = rem(j,nz2); 
  if  g~=0
    X2i(j) = X2ii(1,g);
    Y2i(j) = Y2ii(1,g);
  else
    X2i(j) = X2ii(1,nz2);
    Y2i(j) = Y2ii(1,nz2);
  end
end
X2=X2i.';
Y2=Y2i.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Z1i coordinate   %%%%%%%%%

n=1;
z=-h/2;
for i=1:nz1*hz
    if i <= n*nz1
        n=n;
        z=z;
    else
        n=n+1;
        z=z+Ziv;
    end
    Z1i(i)=z;  
end
%%%%%%%     Z2i coordinate   %%%%%%%%%
n=1;
z=-h/2;
for i=1:nz2*hz
    if i <= n*nz2
        n=n;
        z=z;
    else
        n=n+1;
        z=z+Ziv;
    end
    Z2i(i)=z;  
end
at = 1:1:ni;        %%%atom ID
ab = 1:1:nb;        %%%atom bond ID
b1 = ones(1,ni);
b2 = ones(1,nb);
A1 = at.';
A2 = ab.';
B1 = b1.';
B2 = b2.';
%%%%%% coordinate %%%%%%
X=[X1;X2];
Y=[Y1;Y2];
Z=[Z1i.';Z2i.'];
coordinate = [A1,B1,X,Y,Z];

%%%%%%%%  Bonding  %%%%%%%%
%%%%%%%  R1 crd Bond %%%%%%%%
bo=2:1:ni;
k1=1;
for i=1:1:nz1*hz
    t = rem(i,nz1);
    if t ~=0
        t=t;
        k1=k1;
        Bo(i)=bo(1,i);     %%%Z1為直列bond
    else
        k1=k1+1;
        Bo(i)=1+(k1-2)*nz1;
    end  
end    


%%%%%%% R2 crd Bond %%%%%%%%
bi=nz1*hz+2:1:ni;
k1=1;
for i=1:1:ni-nz1*hz
    t = rem(i,nz2);
    if t ~=0
        t=t;
        k1=k1;
        Bi(i)=bi(1,i);     %%%Z1為直列bond
    else
        k1=k1+1;
        Bi(i)=1+(k1-2+hz)*nz2;
    end  
end    
BD2 = [Bo,Bi];
BD3 = [A1,BD2.'];

%%%%%%%%垂直bond%%%%%%
t3 = 1:1:nz1*(hz-1);
t4 = nz1*hz+1:1:nz*hz-nz2;  %%%層數有在增加的話這裡要修改
t5 = nz1+1:1:nz1*hz;
t6 = nz1*hz+nz2+1:1:nz*hz;
bd4 = [t3.',t5.'];
bd5 = [t4.',t6.'];
BD4 = [bd4;bd5];
BD5 = [BD3;BD4];

%%%%%%%%%%% slash bond %%%%%%%%%%
%%%%%%%%%%% R1層 %%%%%%%%%%%%
t7=2:1:ni-nz;
k2=1;
for i=1:1:nz1*hz-nz1;
    t = rem(i,nz1);
   if t ~=0
        t=t;
        k2=k2;
        bs1(i)=t7(1,i)+nz1;
    else
        k2=k2+1;
        bs1(i)=1+(k2-1)*nz1;       
    end
end
%%%%%%%%%%% R2層 %%%%%%%%%%%%
t8=nz1*hz+2+nz2:1:ni;
k3=1;
for i=1:1:nz2*hz-nz2;
    t = rem(i,nz2);
   if t ~=0
        t=t;
        k3=k3;
        bs2(i)=t8(1,i);
    else
        k3=k3+1;
        bs2(i)=(k3-1)*nz2+nz1*hz+1;       
    end
end





A3=1:1:nz1*hz-nz1;
A4=nz1*hz+1:1:ni-nz2;
A5=[A3,A4];
BS1=bs1.';
BS2=bs2.';
BS=[BS1;BS2];
BS2=[A5.',BS];

%%%%%%%%%%%%%%%%%%%%%

BD=[BD5;BS2];
bond=[A2,B2,BD];

%%%%%%%%   輸出.data file  %%%%%%%%%%%
%%%%%%%%    Header Lines  %%%%%%%%%%%
fileID = fopen('data.0715test','w');

fprintf(fileID,'%1$s\n\n','Lammps Data file');
fprintf(fileID,'%1d\b',ni);
fprintf(fileID,'%1$s\n','atoms');
fprintf(fileID,'%1d\b',nb);
fprintf(fileID,'%1$s\n','bonds');
fprintf(fileID,'%1d\b',nangles);
fprintf(fileID,'%1$s\n','angles');
fprintf(fileID,'%1d\b',ndihedrals);
fprintf(fileID,'%1$s\n','dihedrals');
fprintf(fileID,'%1d\b',nimpropers);
fprintf(fileID,'%1$s\n\n','impropers');

fprintf(fileID,'%1d\b',nt);
fprintf(fileID,'%1$s\n','atom types');
fprintf(fileID,'%1d\b',bt);
fprintf(fileID,'%1$s\n','bond types');
fprintf(fileID,'%1d\b',anglet);
fprintf(fileID,'%1$s\n','angle types');
fprintf(fileID,'%1d\b',dihedralt);
fprintf(fileID,'%1$s\n','dihedral types');
fprintf(fileID,'%1d\b',impropert);
fprintf(fileID,'%1$s\n\n','improper types');

fprintf(fileID,'%1.3f\b %1.3f\b',xo,xi);
fprintf(fileID,'\b%1$s %2$s\n',' xlo ','xhi');
fprintf(fileID,'%1.3f\b %1.3f\b',yo,yi);
fprintf(fileID,'\b%1$s %2$s\n',' ylo ','yhi');
fprintf(fileID,'%1.3f\b %1.3f\b',-zi,zi);
fprintf(fileID,'\b%1$s %2$s\n\n',' zlo ','zhi');

fprintf(fileID,'%1$s\n\n','Masses');
fprintf(fileID,'%1d\b',nt);
fprintf(fileID,'%1.4f\n\n',mass);

%%%%%%%%%%    Atoms    %%%%%%%%
fprintf(fileID,'%1$s\n\n','Atoms');
fprintf(fileID,'%-1d %1d %1.5f %1.5f %1.5f \n',coordinate.');

%%%%%%%%%%    Bonds    %%%%%%%%
fprintf(fileID,'%1$s\n\n','Bonds');
fprintf(fileID,'%-1d %-1d %-1d %-1d \n',bond.');

fclose(fileID);
