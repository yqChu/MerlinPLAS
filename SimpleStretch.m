%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 MERLIN2                               %%
%                     Written by: Ke Liu (ke.liu@gatech.edu)              %
% Ref: K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid    %
%      origami - An efficient computational approach.' PRSA.              %
%      K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to  % 
%      capture highly nonlinear behavior of non-rigid origami.'           %
%      Proceedings of IASS Annual Symposium 2016.                         %
%      E. T. Filipov, K. Liu, T. Tachi, M. Schenk, G. H. Paulino (2017).  %
%      'Bar and hinge models for scalable analysis of origami.'  IJSS     %
%      K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear          %
%      structural analysis of origami assemblages using the MERLIN2       %
%      software.' Origami^7.                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =========== MIURA BENDING =========================================== %%
clear all; close all; clc;
%% Define geomtry
Kf = 0.33; Kb = 3.3; Abar = 2e-1;
% Elastic Modulus
EMod = 10000;
% Yield force & plastic modulus
YBar = 80; Ydf = 100000; Ydb = 100000; 
PMod = 10000; pl_mod_fold = 0.02; pl_mod_bend = 0.2;

Node = [0, 0, 0;
        0, 2, 0;
        1, 1, 0;
        -1, 1, 0;];
Panel = {[1, 2, 3];
         [1, 2, 4];};

% Visualize initial configuration 
figure()
PlotOri(Node,Panel,[],'PanelColor','g')
axis equal; axis off;
light
% Inspect nodal index assignment
figure()
PlotOri(Node,Panel,[],'ShowNumber','on');
axis equal

%% Set up boundary conditions
Supp = [1, 1, 1, 1;
        2, 1, 1, 1;
        4, 1, 1, 1];
Load = [3, 1, 0, 0;];
indp = 3;

%% Define material and modeling parameters
DIRC = [1, 0, 0];
LoadSize = 0.0001;
CycleIcrm = 200;
TotalIcr = 2000;
AnalyInputOpt = struct(...
    'ModelType','N4B5',...
    'MaterCalib','manual',... 
    'MaxIcr', TotalIcr, ...
    'Abar', Abar,...
    'Kb',Kb,...
    'Kf',Kf,...
    'RotSprBend', @(he,h0,Kb,L0)EnhancedLinear(he,h0,Kb,L0,0.1,359.9),...
    'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,0.1,359.9),...
    'LoadType','Displacement',...    % Displacement load
    'AdaptiveLoad', @(Node,U,icrm)LoadFunSimpleStretch(Node,U,icrm,LoadSize,DIRC,CycleIcrm), ...
    'DispStep', TotalIcr);

%% Perform analysis
% Assemble input data
[truss, angles, AnalyInputOpt] = PrepareData(Node,Panel,Supp,Load,AnalyInputOpt,EMod,YBar,PMod,Ydf,Ydb,pl_mod_fold,pl_mod_bend);
% Specify initial deformation state
truss.U0 = zeros(3*size(truss.Node,1),1);
% Perform path-following analysis
[Uhis,Fhis] = PathAnalysis(truss,angles,AnalyInputOpt);
% Postprocess output data
Uhis = real(Uhis);
Fhis = real(Fhis);
STAT = PostProcess(Uhis,truss,angles); 

%% Visualize simulation
instdof = [indp(1),1];
interv = 1; endicrm = size(Uhis,2);
% Animation monitoring node-wise change
VIntensityDataInten = zeros(size(truss.Node,1),size(Uhis,2));
IntensityDataM = bsxfun(@times,STAT.bar.Sx,truss.A);
for k = 1:size(Uhis,2)
    IntensityDataIntenk = sparse(truss.Bars(:,1),truss.Bars(:,2),abs(IntensityDataM(:,k)),size(truss.Node,1),size(truss.Node,1));
    VIntensityDataInten(:,k) = sum((IntensityDataIntenk+IntensityDataIntenk'),2); 
end
% VisualFold(Uhis(:,1:interv:endicrm),truss,angles,-Fhis(1:interv:endicrm,:),instdof,'IntensityMap','Vertex','IntensityData',VIntensityDataInten)
hold off
plot(Uhis(7,1:interv:endicrm),Fhis(1:interv:endicrm,:))
grid on

% VisualFold(Uhis(:,1:interv:endicrm),truss,angles,'none','miura5x5ptd',0.0001,LF_his,instdof,[-inf inf -inf inf])
% Animation monitoring panel-wise value change
% VisualFold(Uhis(:,1:10:endicrm),truss,angles,[],[],'IntensityMap','Edge','IntensityData',STAT.bar.Sx(:,1:10:endicrm),'ShowInitial','off')
% Animation only showing the configurational change
% VisualFold(Uhis(:,1:interv:endicrm),truss,angles,[],[])

%% Plot diagrams
% Load vs displacement
% dsp = sign(instdof(2))*Uhis((instdof(1)*3-(3-abs(instdof(2)))),:);
% figure()
% plot(dsp,Fhis,'b-','linewidth',1);
% axis tight
% xlabel('Displacement','fontsize',14);
% ylabel('Load','fontsize',14);

% Stored energy vs displacement
% figure()
% plot(dsp,STAT.PE,'r-','linewidth',2);    % Red line is the total energy.
% hold on                                  % Between red and cyan is the folding energy. 
% plot(dsp,STAT.bend.UB+STAT.bar.US,'c-'); % Between cyan and magenta is the portion of energy for bending.
% plot(dsp,STAT.bar.US,'m-');              % Below magenta is the stretching energy of bars.
% axis tight
% xlabel('Displacement','fontsize',14);
% ylabel('Stored Energy','fontsize',14);

%% Plot final configuration
% Ux = Uhis(:,end);
% Nodew = truss.Node;
% Nodew(:,1) = truss.Node(:,1)+Ux(1:3:end);
% Nodew(:,2) = truss.Node(:,2)+Ux(2:3:end);
% Nodew(:,3) = truss.Node(:,3)+Ux(3:3:end); 
% figure()
% plot initial configuration
% PlotOri(truss.Node,angles.Panel,truss.Trigl,'FoldEdgeStyle','-','EdgeShade',0.3,'PanelColor','none');
% plot deformed configuration
% PlotOri(Nodew,angles.Panel,truss.Trigl);
% axis equal; axis off;
% camproj('perspective')
% light
% view(117,18)
% rotate3d on

%% Export final configuration to an OBJ file
% Write2OBJ('BendedMiura5x5', Nodew, truss.Trigl, truss.Bars, angles.bend)