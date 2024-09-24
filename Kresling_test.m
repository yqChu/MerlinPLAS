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
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
Kf = 0.2; Kb = 2; Abar = 1e-1;
% Elastic Modulus
EMod = 1e+04;
% Yield force & plastic modulus
% YBar = 800000; Ydf = 100000; Ydb = 1000000; % Elastic 
YBar = 200; Ydf = 0.02; Ydb = 0.2; % Elasto-Plastic
% YBar = 80; Ydf = 0.01; Ydb = 10; % Elasto-Plastic
PMod = 1000; pl_mod_fold = 0.1; pl_mod_bend = 1;

R0 = 1;
theta0 = 11/18*pi;
height0 = 1;
height0 = height0/R0;
Np = 6;
Node = R0*[cos(0), sin(0), 0;
           cos(2*pi/Np), sin(2*pi/Np), 0;
           cos(4*pi/Np), sin(4*pi/Np), 0;
           cos(6*pi/Np), sin(6*pi/Np), 0;
           cos(8*pi/Np), sin(8*pi/Np), 0;
           cos(10*pi/Np), sin(10*pi/Np), 0;
           cos(0+theta0), sin(0+theta0), height0;
           cos(2*pi/Np+theta0), sin(2*pi/Np+theta0), height0;
           cos(4*pi/Np+theta0), sin(4*pi/Np+theta0), height0;
           cos(6*pi/Np+theta0), sin(6*pi/Np+theta0), height0;
           cos(8*pi/Np+theta0), sin(8*pi/Np+theta0), height0;
           cos(10*pi/Np+theta0), sin(10*pi/Np+theta0), height0;];
Panel = {[1, 2, 7];
         [2, 3, 8];
         [3, 4, 9];
         [4, 5, 10];
         [5, 6, 11];
         [6, 1, 12];
         [7, 8, 2];
         [8, 9, 3];
         [9, 10, 4];
         [10, 11, 5];
         [11, 12, 6];
         [12, 7, 1]};

% Visualize initial configuration 
figure()
PlotOri(Node,Panel,[],'PanelColor','g')
axis equal; axis off;
light
% Inspect nodal index assignment
figure()
PlotOri(Node,Panel,[],'ShowNumber','on');
axis equal; axis off;

%% Set up boundary conditions
m = size(Node,1);
Supp = [1, 1, 1, 1;
        2, 1, 1, 1;
        3, 1, 1, 1;
        4, 1, 1, 1;
        5, 1, 1, 1;
        6, 1, 1, 1;];
Load = [ 5, 0, 0, 1;
         6, 0, 0, 1;
         7, 0, 0, 1;
         8, 0, 0, 1];
indp = 7;

%% Define material and modeling parameters

% Simulation options using the N4B5 model
AnalyInputOpt = struct(...
    'ModelType','N4B5',...
    'MaterCalib','manual',... 
    'MaxIcr', 1000, ...
    'Abar', Abar,...
    'Kb',Kb,...
    'Kf',Kf,...
    'RotSprBend', @(he,h0,Kb,L0)EnhancedLinear(he,h0,Kb,L0,30,330),...
    'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,30,330),...
    'LoadType','Displacement',...    % Displacement load
    'AdaptiveLoad',@LoadFunKresling,...
    'DispStep', 1000);

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
instdof = [indp(1),3];
interv = 1; endicrm = size(Uhis,2);
% Animation monitoring node-wise change
VIntensityDataInten = zeros(size(truss.Node,1),size(Uhis,2));
IntensityDataM = bsxfun(@times,STAT.bar.Sx,truss.A);
for k = 1:size(Uhis,2)
    IntensityDataIntenk = sparse(truss.Bars(:,1),truss.Bars(:,2),abs(IntensityDataM(:,k)),size(truss.Node,1),size(truss.Node,1));
    VIntensityDataInten(:,k) = sum((IntensityDataIntenk+IntensityDataIntenk'),2); 
end
VisualFold(Uhis(:,1:interv:endicrm),truss,angles,-Fhis(1:interv:endicrm,:),instdof,'IntensityMap','Vertex','IntensityData',VIntensityDataInten)
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
Ux = Uhis(:,end);
Nodew = truss.Node;
Nodew(:,1) = truss.Node(:,1)+Ux(1:3:end);
Nodew(:,2) = truss.Node(:,2)+Ux(2:3:end);
Nodew(:,3) = truss.Node(:,3)+Ux(3:3:end); 
figure()
hold off
% plot initial configuration
PlotOri(truss.Node,angles.Panel,truss.Trigl,'FoldEdgeStyle','-','EdgeShade',0.3,'PanelColor','none');
hold off
% plot deformed configuration
PlotOri(Nodew,angles.Panel,truss.Trigl);
axis equal; axis off;
camproj('perspective')
light
view(117,18)
rotate3d on

%% Export final configuration to an OBJ file
% Write2OBJ('BendedMiura5x5', Nodew, truss.Trigl, truss.Bars, angles.bend)