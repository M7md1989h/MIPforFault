  clear all
      clc
load P
%    pwithDG to get the measured voltage index
k = find(P);
%% Y matrix
% Initialize OpenDSS
% Create the OpenDSS Object
DSSobj = actxserver('OpenDSSEngine.DSS');
% Start up the solver
if ~DSSobj.Start(0),
disp('Unable to start the OpenDSS Engine')
return
end
% Set up the interface variables
DSSText = DSSobj.Text;
DSSCircuit=DSSobj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;
% Run your OpenDSS file 
DSSText.command='Compile (134.dss)'; 
 DSSText.command='vsource.source.enabled=no';
 DSSText.command='Generator.pv134.enabled=no';
 DSSText.command='Generator.pv99.enabled=no';
%  DSSText.command=' batchedit load..* enabled=noY';
DSSSolution.Solve;
% Y matrix
n = length(DSSCircuit.AllNodeNames);
ineven=2:2:n*2;
inodd=1:2:(n*2-1);
Ybus=DSSCircuit.SystemY;
Y=reshape(Ybus,n*2,n)';
Y=Y(:,inodd)+1i*Y(:,ineven);
 YY=Y(4:n-6,4:end);

%
ab=zeros(3,n-3);
ab(1,n-5-3)=1;
ab(2,n-4-3)=1;
ab(3,n-3-3)=1;
ac=zeros(3,n-3);
ac(1,n-2-3)=1;
ac(2,n-1-3)=1;
ac(3,n-3)=1;

 YY=[YY;ab;ac];
ZZ=inv(YY);
%% this section is used to find the measured pre-Voltages(on 15 buses)
% Initialize OpenDSS
% Create the OpenDSS Object
DSSobj = actxserver('OpenDSSEngine.DSS');
% Start up the solver
if ~DSSobj.Start(0),
disp('Unable to start the OpenDSS Engine')
return
end
% Set up the interface variables
DSSText = DSSobj.Text;
DSSCircuit=DSSobj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;
% Run your OpenDSS file 
DSSText.command='Compile (134.dss)'; 
% DSSText.command='solve';
% DSSSolution.Solve;
n = length(DSSCircuit.AllNodeNames);
ineven=2:2:n*2;
inodd=1:2:(n*2-1);
% pre-fault Voltage measurments
vckt= DSSCircuit.YNodeVarray';
VVbefore =vckt(inodd)+1i*vckt(ineven);

% %% this section is used to find the measured post-Voltages(on 15 buses)
% % % % 
% % 
for s=2:134
DSSobj = actxserver('OpenDSSEngine.DSS');
% Start up the solver
if ~DSSobj.Start(0),
disp('Unable to start the OpenDSS Engine')
return
end
% Set up the interface variables
DSSText = DSSobj.Text;
DSSCircuit=DSSobj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;
% Run your OpenDSS file 
DSSText.command='Compile (134.dss)';%place the OpenDss file path here 
DSSText.command='Solve mode=dynamics stepsize=0.00002';

formatSpec='new fault.L phases=1 bus1 =%s r=0.5';

 A1=DSSCircuit.AllBusNames{s+1}

 fault=sprintf(formatSpec,A1);
  DSSText.command=fault;
DSSText.command='solve';
DSSSolution.Solve;
nodes=DSSCircuit.YNodeOrder';
% Voltage
vckt= DSSCircuit.YNodeVarray';
VVduring =vckt(inodd)+1i*vckt(ineven);

%% after collecting the Data, we use this section to analysis and find the faulted bus. 
Vbefore=VVbefore(4:end);
Vduring=VVduring(4:end);
deltaVmeasured=(Vduring-Vbefore);
M=1e6;
 
cvx_solver mosek
cvx_begin quiet
 
%                         cvx_precision best

variable Id(n-9) complex
variable u(n-9) binary
variable I_F(n-9) complex
expression vv(n-3) 
 variable y(134) binary


vv=ZZ*[-Id;deltaVmeasured(n-5-3);deltaVmeasured(n-4-3);deltaVmeasured(n-3-3);deltaVmeasured(n-2-3);deltaVmeasured(n-1-3);deltaVmeasured(n-3)];

subject to
   sum(y)==1
for i=1:n-9
      -(1-u(i))*M + real(I_F(i)) <= real(Id(i)) <= real(I_F(i)) + (1-u(i))*M
      -(1-u(i))*M + imag(I_F(i)) <= imag(Id(i)) <= imag(I_F(i)) + (1-u(i))*M
         -u(i)*M <= real(Id(i)) <=  u(i)*M
         -u(i)*M  <= imag(Id(i)) <=  u(i)*M
end


c=[1:134];
C=[1:3:n-9];
u(C)+u(C+1)+u(C+2)==3*(y(c))
        y(1)==0
      

expression cost
cost =norm(deltaVmeasured(k)-vv(k));  


minimize cost

cvx_end 


fault_bus = find(u>0);
% faultedbus=nodes(fault_bus)
faultedbus=nodes(fault_bus+3);
Ifault=abs(Id(fault_bus));
Ifaultangle=angle(Id(fault_bus));
%    degree=Ifaultangle/pi*180
faultedphase0(s,1)={A1};
faultedphase0(s,2)=faultedbus(1,1);
faultedphase0(s,3)=faultedbus(1,2);
faultedphase0(s,4)=faultedbus(1,3);
faultedphase0(s,5)={Ifault(1,1)};
faultedphase0(s,6)={Ifault(2,1)};
faultedphase0(s,7)={Ifault(3,1)};
faultedphase0(s,8)={cvx_optval};
end
%  save faultedphase0
% 
 for s=2:134
DSSobj = actxserver('OpenDSSEngine.DSS');
% Start up the solver
if ~DSSobj.Start(0),
disp('Unable to start the OpenDSS Engine')
return
end
% Set up the interface variables
DSSText = DSSobj.Text;
DSSCircuit=DSSobj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;
% Run your OpenDSS file 
DSSText.command='Compile (134.dss)';%place the OpenDss file path here 
DSSText.command='Solve mode=dynamics stepsize=0.00002';

formatSpec='new fault.LL phases=2 bus1 =%s r=0.5';

 A1=DSSCircuit.AllBusNames{s+1}

 fault=sprintf(formatSpec,A1);
  DSSText.command=fault;
DSSText.command='solve';
DSSSolution.Solve;
nodes=DSSCircuit.YNodeOrder';
% Voltage
vckt= DSSCircuit.YNodeVarray';
VVduring =vckt(inodd)+1i*vckt(ineven);

%% after collecting the Data, we use this section to analysis and find the faulted bus. 
Vbefore=VVbefore(4:end);
Vduring=VVduring(4:end);
deltaVmeasured=(Vduring-Vbefore);
M=1e6;
 
cvx_solver mosek
cvx_begin quiet
 
%                         cvx_precision best

variable Id(n-9) complex
variable u(n-9) binary
variable I_F(n-9) complex
expression vv(n-3) 
 variable y(134) binary

vv=ZZ*[-Id;deltaVmeasured(n-5-3);deltaVmeasured(n-4-3);deltaVmeasured(n-3-3);deltaVmeasured(n-2-3);deltaVmeasured(n-1-3);deltaVmeasured(n-3)];

subject to
   sum(y)==1
for i=1:n-9
      -(1-u(i))*M + real(I_F(i)) <= real(Id(i)) <= real(I_F(i)) + (1-u(i))*M
      -(1-u(i))*M + imag(I_F(i)) <= imag(Id(i)) <= imag(I_F(i)) + (1-u(i))*M
         -u(i)*M <= real(Id(i)) <=  u(i)*M
         -u(i)*M  <= imag(Id(i)) <=  u(i)*M
end


c=[1:134];
C=[1:3:n-9];
u(C)+u(C+1)+u(C+2)==3*(y(c))
        y(1)==0
      

expression cost
cost =norm(deltaVmeasured(k)-vv(k));  

minimize cost

cvx_end 


fault_bus = find(u>0);
% faultedbus=nodes(fault_bus)
faultedbus=nodes(fault_bus+3);
Ifault=abs(Id(fault_bus));
Ifaultangle=angle(Id(fault_bus));
%    degree=Ifaultangle/pi*180
faultedphase2(s,1)={A1};
faultedphase2(s,2)=faultedbus(1,1);
faultedphase2(s,3)=faultedbus(1,2);
faultedphase2(s,4)=faultedbus(1,3);
faultedphase2(s,5)={Ifault(1,1)};
faultedphase2(s,6)={Ifault(2,1)};
faultedphase2(s,7)={Ifault(3,1)};
faultedphase2(s,8)={cvx_optval};
faultedphase2;
 end
%  
%  save faultedphase2

 for s=2:134
DSSobj = actxserver('OpenDSSEngine.DSS');
% Start up the solver
if ~DSSobj.Start(0),
disp('Unable to start the OpenDSS Engine')
return
end
% Set up the interface variables
DSSText = DSSobj.Text;
DSSCircuit=DSSobj.ActiveCircuit;
DSSSolution=DSSCircuit.Solution;
% Run your OpenDSS file 
DSSText.command='Compile (134.dss)';%place the OpenDss file path here 
DSSText.command='Solve mode=dynamics stepsize=0.00002';

formatSpec='new fault.LLL phases=3 bus1 =%s r=0.5';

 A1=DSSCircuit.AllBusNames{s+1}

 fault=sprintf(formatSpec,A1);
  DSSText.command=fault;
DSSText.command='solve';
DSSSolution.Solve;
nodes=DSSCircuit.YNodeOrder';
% Voltage
vckt= DSSCircuit.YNodeVarray';
VVduring =vckt(inodd)+1i*vckt(ineven);

%% after collecting the Data, we use this section to analysis and find the faulted bus. 
Vbefore=VVbefore(4:end);
Vduring=VVduring(4:end);
deltaVmeasured=(Vduring-Vbefore);
M=1e6;
 
cvx_solver mosek
cvx_begin quiet
 
%                         cvx_precision best

variable Id(n-9) complex
variable u(n-9) binary
variable I_F(n-9) complex
expression vv(n-3) 
 variable y(134) binary


vv=ZZ*[-Id;deltaVmeasured(n-5-3);deltaVmeasured(n-4-3);deltaVmeasured(n-3-3);deltaVmeasured(n-2-3);deltaVmeasured(n-1-3);deltaVmeasured(n-3)];

subject to
   sum(y)==1
for i=1:n-9
      -(1-u(i))*M + real(I_F(i)) <= real(Id(i)) <= real(I_F(i)) + (1-u(i))*M
      -(1-u(i))*M + imag(I_F(i)) <= imag(Id(i)) <= imag(I_F(i)) + (1-u(i))*M
         -u(i)*M <= real(Id(i)) <=  u(i)*M
         -u(i)*M  <= imag(Id(i)) <=  u(i)*M
end


c=[1:134];
C=[1:3:n-9];
u(C)+u(C+1)+u(C+2)==3*(y(c))
        y(1)==0
      

expression cost
%  cost =norm(VVduring(k)-vv(k));  
cost =norm(deltaVmeasured(k)-vv(k));  


minimize cost

cvx_end 


fault_bus = find(u>0);
% faultedbus=nodes(fault_bus)
faultedbus=nodes(fault_bus+3);
Ifault=abs(Id(fault_bus));
Ifaultangle=angle(Id(fault_bus));
%    degree=Ifaultangle/pi*180
faultedphase3(s,1)={A1};
faultedphase3(s,2)=faultedbus(1,1);
faultedphase3(s,3)=faultedbus(1,2);
faultedphase3(s,4)=faultedbus(1,3);
faultedphase3(s,5)={Ifault(1,1)};
faultedphase3(s,6)={Ifault(2,1)};
faultedphase3(s,7)={Ifault(3,1)};
faultedphase3(s,8)={cvx_optval};
faultedphase3;
 end
%  save faultedphase3

% S1=find((str2double(faultedphase(:,1 ))-(str2double(faultedphase(:,2 ))) < 0.4));
% S2=find((str2double(faultedphase1(:,1 ))-(str2double(faultedphase1(:,2 )))<  0.4));
% S3=find((str2double(faultedphase2(:,1 ))-(str2double(faultedphase2(:,2 )))<  0.4));
% S4=find((str2double(faultedphase10(:,1 ))-(str2double(faultedphase10(:,2 )))<  0.4));
% S5=find((str2double(faultedphase11(:,1 ))-(str2double(faultedphase11(:,2 )))<  0.4));
% S6=find((str2double(faultedphase12(:,1 ))-(str2double(faultedphase12(:,2 )))<  0.4));
