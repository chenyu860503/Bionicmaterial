%%%%%%%%   lammps data file  %%%%%%%%%%
clc
clear
close all

%%%%%%%%   define parameters %%%%%%%%%%
mass=28;     %%原子量
ri=8;   %%%內層r
ro=10;  %%%外層r
l=2;    %%l為層數
nt=1;   %%%nt為atom種類
bt=1;   %%%bt為bond種類;  %%ni=atoms數目 
nz=8;   %%nz為單位高lattice數目
hz=10;  %%hz為垂直高度原子數目
Ziv=3;     %%上下原子間隔差
h=Ziv*hz;    %%圓柱高度
ni=l*nz*hz; %%atoms數目
theta1=0;
theta2=2*pi/nz;
theta3=theta1:theta2:2*pi-theta2;
nb=ni+(2+l)*(hz-1)*nz; %%nb為一層總bond數目 
nangles=0;      %%angles數目
ndihedrals=0;   %%dihedrals數目
nimpropers=0;   %%impropers數目
anglet=0;     %%angles type
dihedralt=0;  %%dihedral type
impropert=0;  %%improper type



%%%%%%%%%%   box 大小 %%%%%%%%
xo=-200;
xi=200;
yo=-200;
yi=200;
zo=-200;
zi=200;

%%%%%%%atom cooridnate%%%%%%%
%%%%%%%     Xi,Yi建立   %%%%%%%
Xii = [ri*cos(theta3)];
Yii = [ri*sin(theta3)];  
at = 1:1:ni;        %%%atom ID
ab = 1:1:nb;        %%%atom bond ID
n = 1;
for j=1:1:ni/l
  g = rem(j,nz); 
  if  g~=0
    X1(j) = Xii(1,g);
    Y1(j) = Yii(1,g);
  else
    X1(j) = Xii(1,nz);
    Y1(j) = Yii(1,nz);
  end
end
Xi=X1.';
Yi=Y1.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Zi coordinate   %%%%%%%%%
b1 = ones(1,ni);
b2 = ones(1,nb);
A1 = at.';
A2 = ab.';
B1 = b1.';
B2 = b2.';
z=1;
for i=1:ni/l
    if i <= n*nz
        n=n;
        z=z;
    else
        n=n+1;
        z=z+Ziv;
    end
    Zi(i)=z;  
end


%%%%%%%  外層crd Bond %%%%%%%%
bo=2:1:ni;
k1=1;
for i=1:1:ni/l
    t = rem(i,nz);
    if t ~=0
        t=t;
        k1=k1;
        Bo(i)=bo(1,i);     %%%Z1為直列bond
    else
        k1=k1+1;
        Bo(i)=1+(k1-2)*nz;
    end  
end    

%%%%%%%     Xo,Yo建立   %%%%%%%
Xoo = [ro*cos(theta3)];
Yoo = [ro*sin(theta3)];  
n = 1;
for j=1:1:ni/l
  g = rem(j,nz); 
  if  g~=0
    X2(j) = Xoo(1,g);
    Y2(j) = Yoo(1,g);
  else
    X2(j) = Xoo(1,nz);
    Y2(j) = Yoo(1,nz);
  end
end
Xo=X2.';
Yo=Y2.';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     Zo coordinate   %%%%%%%%%
zo=1;
for i=1:ni/l
    if i <= n*nz
        n=n;
        zo=zo;
    else
        n=n+1;
        zo=zo+Ziv;
    end
    Zo(i)=zo;  
end

%%%%%% coordinate %%%%%%
X=[Xi;Xo].*0.001;
Y=[Yi;Yo].*0.001;
Z=[Zi.';Zo.'].*0.001;
coordinate = [A1,B1,X,Y,Z];

%%%%%%%%  Bonding  %%%%%%%%

%%%%%%% 內層crd Bond %%%%%%%%
bi=nz*hz+2:1:ni;
k1=1;
for i=1:1:ni/l
    t = rem(i,nz);
    if t ~=0
        t=t;
        k1=k1;
        Bi(i)=bi(1,i);     %%%Z1為直列bond
    else
        k1=k1+1;
        Bi(i)=1+(k1-2+hz)*nz;
    end  
end    
BD2 = [Bo,Bi];
BD3 = [A1,BD2.'];

%%%%%%%%垂直bond%%%%%%
t3 = 1:1:nz*(hz-1);
t4 = nz*hz+1:1:nz*(l*hz-1);  %%%層數有在增加的話這裡要修改
t5 = nz+1:1:nz*hz;
t6 = nz*(hz+1)+1:1:nz*l*hz;
bd4 = [t3.',t5.'];
bd5 = [t4.',t6.'];
BD4 = [bd4;bd5];
BD5 = [BD3;BD4];

%%%%%% slash 內層ˋ外層bond %%%%%%
t7=2:1:ni-nz;
k2=1;
for i=1:1:ni/l-nz;
    t = rem(i,nz)
    if t ~=0
        t=t;
        k2=k2;
        bso(i)=t7(1,i)+nz;
        bsi(i)=ni/l+t7(1,i)+nz;
    else
        k2=k2+1;
        bso(i)=1+(k2-1)*nz;
        bsi(i)=ni/l+1+(k2-1)*nz;
    end
end
A3=1:1:ni/l-nz;
A4=ni/l+1:1:ni-nz;
A5=[A3,A4];

BSo=bso.';
BSi=bsi.';
BS=[BSo;BSi];
BS2=[A5.',BS];

%%%%%%%%%%%%%%%%%%%%%

BD=[BD5;BS2];
bond=[A2,B2,BD];

%%%%%%%%   輸出.data file  %%%%%%%%%%%
%%%%%%%%    Header Lines  %%%%%%%%%%%
fileID = fopen('data.cyldslashMultiLayer100','w');

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
fprintf(fileID,'%1.3f\b %1.3f\b',zo,zi);
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
