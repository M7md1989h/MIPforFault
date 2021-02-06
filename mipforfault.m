 clear all
  clc
load P %to get the measured voltage index
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
DSSText.command='Compile (Run_IEEE123Bus.DSS)'; 
 DSSText.command='vsource.source.enabled=no';
%DSSText.command=' batchedit load..* enabled=noY'; % Load data are not included in the Ybus if this line is enable.  
DSSText.command='solve';
DSSSolution.Solve;
% Y matrix
n = length(DSSCircuit.AllNodeNames);

C=str2double(DSSCircuit.AllBusNames(2:end));
CC=str2double(DSSCircuit.YNodeOrder(4:end));

ineven=2:2:n*2;
inodd=1:2:(n*2-1);
Ybus=DSSCircuit.SystemY;
Y=reshape(Ybus,n*2,n)';
Y=Y(:,inodd)+1i*Y(:,ineven);
 Z=inv(Y);

YY=Y(4:end,4:end);
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
DSSText.command='Compile (Run_IEEE123Bus.DSS)'; 
DSSText.command='solve';
DSSSolution.Solve;
n = length(DSSCircuit.AllNodeNames);
ineven=2:2:n*2;
inodd=1:2:(n*2-1);
% pre-fault Voltage measurments
vckt= DSSCircuit.YNodeVarray';
VVbefore =vckt(inodd)+1i*vckt(ineven);

%%
for i=1:129
idx = find( abs(CC(:,1 ) - C(i)) <  0.4);
L=length(idx);
if L==3
X(i)=C(i);

elseif L==2
  XX(i)=C(i);
else 
    XXX(i)=C(i);
end
end

threephasebus=nonzeros(X);
twophasebus=nonzeros(XX);
onephasebus=nonzeros(XXX);
three_twophasebus=[threephasebus;twophasebus];
%%
  for s=1:269
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
DSSText.command='Compile (Run_IEEE123Bus.DSS)';%place the OpenDss file path here 
formatSpec='new fault.L phases=1 bus1 =%s r=0';
A1=DSSCircuit.AllNodeNames{s+3}
% % formatSpec='new fault.LL phases=2 bus1 = %d.2.3';
%   A1=three_twophasebus(s)
%   A2=three_twophasebus(s);
% % formatSpec='new fault.LLL phases=3 bus1 = %d';
% A1=threephasebus(s)
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
deltaVmeasured=Vduring-Vbefore;

% %%noise
% % Vbeforenoise(k,1)=((abs(Vbefore(k))).*(1+(random('norm', 0, 0.005, size(k) )))).*exp(j.*(angle(Vbefore(k))));
% % Vduringnoise(k,1)=((abs(Vduring(k))).*(1+(random('norm', 0, 0.005, size(k) )))).*exp(j.*(angle(Vduring(k))));
% % deltaVmeasured=Vduringnoise-Vbeforenoise;


M=1e6;
  
% % 
cvx_solver mosek
 
cvx_begin quiet
 
%                           cvx_precision high

variable Id(n-3) complex
variable u(n-3) binary
variable I_F(n-3) complex
expression vv(n-3) 
       variable y(129) binary




vv=ZZ*(-Id);

subject to
%                sum(u)==1

      sum(y)==1
for i=1:n-3
      -(1-u(i))*M + real(I_F(i)) <= real(Id(i)) <= real(I_F(i)) + (1-u(i))*M
      -(1-u(i))*M + imag(I_F(i)) <= imag(Id(i)) <= imag(I_F(i)) + (1-u(i))*M
         -u(i)*M <= real(Id(i)) <=  u(i)*M
         -u(i)*M  <=imag(Id(i)) <=  u(i)*M
end
% 
% 
% 
for i=1:129
idx = find( abs(CC(:,1 ) - C(i)) <  0.4);
L=length(idx);
if L==3
 u(idx(1))+u(idx(2))+u(idx(3))==3*(y(i))

elseif L==2
  
 u(idx(1))+u(idx(2))==2*(y(i))



    
else 
 u(idx(1))==1*(y(i))

end
end
       
expression cost
%  cost =norm(VVduring(k)-vv(k));  
cost =norm(deltaVmeasured(k)-vv(k));  


minimize cost

cvx_end 

fault_bus = find(u==1);
% faultedbus=nodes(fault_bus)
faultedbus=nodes(fault_bus+3);
Ifault=abs(Id(fault_bus));
Ifaultangle=angle(Id(fault_bus));

if length(fault_bus)==3
faultedphase(s,1)={A1};
faultedphase(s,2)=faultedbus(1,1);
faultedphase(s,3)=faultedbus(1,2);
faultedphase(s,4)=faultedbus(1,3);
faultedphase(s,5)={Ifault(1,1)};
faultedphase(s,6)={Ifault(2,1)};
faultedphase(s,7)={Ifault(3,1)};
%  faultedphase(s,9)={(Ifaultangle(1,1))*180/pi};
%  faultedphase(s,10)={(Ifaultangle(2,1))*180/pi};
%  faultedphase(s,11)={(Ifaultangle(3,1))*180/pi};

faultedphase(s,8)={cvx_optval};
elseif length(fault_bus)==2
    faultedphase(s,1)={A1};
faultedphase(s,2)=faultedbus(1,1);
faultedphase(s,3)=faultedbus(1,2);
faultedphase(s,4)=faultedbus(1,2);
faultedphase(s,5)={Ifault(1,1)};
faultedphase(s,6)={Ifault(2,1)};
faultedphase(s,7)={Ifault(2,1)};
faultedphase(s,8)={cvx_optval};
% faultedphase(s,9)={(Ifaultangle(1,1))*180/pi};
% faultedphase(s,10)={(Ifaultangle(2,1))*180/pi};
% faultedphase(s,11)={(Ifaultangle(2,1))*180/pi};
else 
faultedphase(s,1)={A1};
faultedphase(s,2)=faultedbus(1,1);
faultedphase(s,3)=faultedbus(1,1);
faultedphase(s,4)=faultedbus(1,1);
faultedphase(s,5)={Ifault(1,1)};
faultedphase(s,6)={Ifault(1,1)};
faultedphase(s,7)={Ifault(1,1)};
faultedphase(s,8)={cvx_optval};
% faultedphase(s,9)={(Ifaultangle(1,1))*180/pi};
% faultedphase(s,10)={(Ifaultangle(1,1))*180/pi};
% faultedphase(s,11)={(Ifaultangle(1,1))*180/pi};
end
  end

% % % % % % % %   
  for s=1:72
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
DSSText.command='Compile (Run_IEEE123Bus.DSS)';%place the OpenDss file path here 

% formatSpec='new fault.L phases=1 bus1 =%s';
% A1=DSSCircuit.AllNodeNames{s+3}
formatSpec='new fault.LL phases=2 bus1 = %d r=0';
A1=three_twophasebus(s)
% % formatSpec='new fault.LLL phases=3 bus1 = %d';
% % A1=threephasebus(s)
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
    deltaVmeasured=Vduring-Vbefore;
% % % %%noise
% Vbeforenoise(k,1)=((abs(Vbefore(k))).*(1+(random('norm', 0, 0.005, size(k) )))).*exp(j.*(angle(Vbefore(k))));
% Vduringnoise(k,1)=((abs(Vduring(k))).*(1+(random('norm', 0, 0.005, size(k) )))).*exp(j.*(angle(Vduring(k))));
% deltaVmeasured=Vduringnoise-Vbeforenoise;

M=1e6;
  
% % 
cvx_solver mosek
 
cvx_begin quiet
 
%         cvx_precision high

variable Id(n-3) complex
variable u(n-3) binary
variable I_F(n-3) complex
expression vv(n-3) 
variable y(129) binary




vv=ZZ*(-Id);

subject to

 sum(y)==1
for i=1:n-3
      -(1-u(i))*M + real(I_F(i)) <= real(Id(i)) <= real(I_F(i)) + (1-u(i))*M
      -(1-u(i))*M + imag(I_F(i)) <= imag(Id(i)) <= imag(I_F(i)) + (1-u(i))*M
         -u(i)*M <= real(Id(i)) <=  u(i)*M
         -u(i)*M  <=imag(Id(i)) <=  u(i)*M
end



for i=1:129
idx = find( abs(CC(:,1 ) - C(i)) <  0.4);
L=length(idx);
if L==3
 u(idx(1))+u(idx(2))+u(idx(3))==3*(y(i))

elseif L==2
  
 u(idx(1))+u(idx(2))==2*(y(i))



    
else 
 u(idx(1))==1*(y(i))

end
end
expression cost
%  cost =norm(VVduring(k)-vv(k));  
cost =norm(deltaVmeasured(k)-vv(k));  


minimize cost

cvx_end 


fault_bus = find(u==1);
% faultedbus=nodes(fault_bus)
faultedbus=nodes(fault_bus+3);
Ifault=abs(Id(fault_bus));
Ifaultangle=angle(Id(fault_bus));

if length(fault_bus)==3
faultedphase0(s,1)={A1};
faultedphase0(s,2)=faultedbus(1,1);
faultedphase0(s,3)=faultedbus(1,2);
faultedphase0(s,4)=faultedbus(1,3);
faultedphase0(s,5)={Ifault(1,1)};
faultedphase0(s,6)={Ifault(2,1)};
faultedphase0(s,7)={Ifault(3,1)};
faultedphase0(s,8)={cvx_optval};
faultedphase0(s,9)={(Ifaultangle(1,1))*180/pi};
 faultedphase0(s,10)={(Ifaultangle(2,1))*180/pi};
 faultedphase0(s,11)={(Ifaultangle(3,1))*180/pi};
elseif length(fault_bus)==2
    faultedphase0(s,1)={A1};
faultedphase0(s,2)=faultedbus(1,1);
faultedphase0(s,3)=faultedbus(1,2);
faultedphase0(s,4)=faultedbus(1,2);
faultedphase0(s,5)={Ifault(1,1)};
faultedphase0(s,6)={Ifault(2,1)};
faultedphase0(s,7)={Ifault(2,1)};
faultedphase0(s,8)={cvx_optval};
faultedphase0(s,9)={(Ifaultangle(1,1))*180/pi};
 faultedphase0(s,10)={(Ifaultangle(2,1))*180/pi};
 faultedphase0(s,11)={(Ifaultangle(2,1))*180/pi};
else 
faultedphase0(s,1)={A1};
faultedphase0(s,2)=faultedbus(1,1);
faultedphase0(s,3)=faultedbus(1,1);
faultedphase0(s,4)=faultedbus(1,1);
faultedphase0(s,5)={Ifault(1,1)};
faultedphase0(s,6)={Ifault(1,1)};
faultedphase0(s,7)={Ifault(1,1)};
faultedphase0(s,8)={cvx_optval};
faultedphase0(s,9)={(Ifaultangle(1,1))*180/pi};
 faultedphase0(s,10)={(Ifaultangle(2,1))*180/pi};
 faultedphase0(s,11)={(Ifaultangle(3,1))*180/pi};
end
  end
%   
  for s=1:68
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
DSSText.command='Compile (Run_IEEE123Bus.DSS)';%place the OpenDss file path here 

% formatSpec='new fault.L phases=1 bus1 =%s';
% A1=DSSCircuit.AllNodeNames{s+3}
% formatSpec='new fault.LL phases=2 bus1 = %d';
%  A1=three_twophasebus(s)
formatSpec='new fault.LLL phases=3 bus1 = %d r=0';
A1=threephasebus(s)
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
     deltaVmeasured=Vduring-Vbefore;

% % %%noise
% Vbeforenoise(k,1)=((abs(Vbefore(k))).*(1+(random('norm', 0, 0.005, size(k) )))).*exp(j.*(angle(Vbefore(k))));
% Vduringnoise(k,1)=((abs(Vduring(k))).*(1+(random('norm', 0, 0.005, size(k) )))).*exp(j.*(angle(Vduring(k))));
% deltaVmeasured=Vduringnoise-Vbeforenoise;



M=1e6;
  
% % 
cvx_solver mosek
 
cvx_begin quiet
 
%         cvx_precision high

variable Id(n-3) complex
variable u(n-3) binary
variable I_F(n-3) complex
expression vv(n-3) 
variable y(129) binary




vv=ZZ*(-Id);

subject to
%       sum(u)==1

 sum(y)==1
for i=1:n-3
      -(1-u(i))*M + real(I_F(i)) <= real(Id(i)) <= real(I_F(i)) + (1-u(i))*M
      -(1-u(i))*M + imag(I_F(i)) <= imag(Id(i)) <= imag(I_F(i)) + (1-u(i))*M
         -u(i)*M <= real(Id(i)) <=  u(i)*M
         -u(i)*M  <=imag(Id(i)) <=  u(i)*M
end



for i=1:129
idx = find( abs(CC(:,1 ) - C(i)) <  0.4);
L=length(idx);
if L==3
 u(idx(1))+u(idx(2))+u(idx(3))==3*(y(i))

elseif L==2
  
 u(idx(1))+u(idx(2))==2*(y(i))



    
else 
 u(idx(1))==1*(y(i))

end
end
expression cost
%  cost =norm(VVduring(k)-vv(k));  
cost =norm(deltaVmeasured(k)-vv(k));  


minimize cost

cvx_end 


fault_bus = find(u==1);
% faultedbus=nodes(fault_bus)
faultedbus=nodes(fault_bus+3);
Ifault=abs(Id(fault_bus));
Ifaultangle=angle(Id(fault_bus));

if length(fault_bus)==3
faultedphase1(s,1)={A1};
faultedphase1(s,2)=faultedbus(1,1);
faultedphase1(s,3)=faultedbus(1,2);
faultedphase1(s,4)=faultedbus(1,3);
faultedphase1(s,5)={Ifault(1,1)};
faultedphase1(s,6)={Ifault(2,1)};
faultedphase1(s,7)={Ifault(3,1)};
faultedphase1(s,8)={cvx_optval};

% faultedphase1(s,9)={(Ifaultangle(1,1))*180/pi};
%  faultedphase1(s,10)={(Ifaultangle(2,1))*180/pi};
%  faultedphase1(s,11)={(Ifaultangle(3,1))*180/pi};

elseif length(fault_bus)==2
    faultedphase1(s,1)={A1};
faultedphase1(s,2)=faultedbus(1,1);
faultedphase1(s,3)=faultedbus(1,2);
faultedphase1(s,4)=faultedbus(1,2);
faultedphase1(s,5)={Ifault(1,1)};
faultedphase1(s,6)={Ifault(2,1)};
faultedphase1(s,7)={Ifault(2,1)};
faultedphase1(s,8)={cvx_optval};
% faultedphase1(s,9)={(Ifaultangle(1,1))*180/pi};
%  faultedphase1(s,10)={(Ifaultangle(2,1))*180/pi};
%  faultedphase1(s,11)={(Ifaultangle(2,1))*180/pi};

else 
faultedphase1(s,1)={A1};
faultedphase1(s,2)=faultedbus(1,1);
faultedphase1(s,3)=faultedbus(1,1);
faultedphase1(s,4)=faultedbus(1,1);
faultedphase1(s,5)={Ifault(1,1)};
faultedphase1(s,6)={Ifault(1,1)};
faultedphase1(s,7)={Ifault(1,1)};
faultedphase1(s,8)={cvx_optval};
% faultedphase1(s,9)={(Ifaultangle(1,1))*180/pi};
%  faultedphase1(s,10)={(Ifaultangle(1,1))*180/pi};
%  faultedphase1(s,11)={(Ifaultangle(1,1))*180/pi};

end
  end
  
% % % %   
%  S=find(str2double(faultedphase(:,1))==str2double(faultedphase(:,2))|str2double(faultedphase(:,1))==str2double(faultedphase(:,3))| ...
%     str2double(faultedphase(:,1))==str2double(faultedphase(:,4)));
