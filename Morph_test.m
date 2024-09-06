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
sec_hor=4;  sec_vert=4; 
a = 2; b = 2; alfa1 = 60; alfa2 = 40.02;
beta = 40; sgngamma2 = 40;
% Material-related parameters: Kf-folding stiffness; Kb-bending stiffness;
Kf = 0.0033; Kb = 0.033; Abar = 2e-5;
% E0-stretching stiffness; Abar-bar area (assume uniform)
limlft = 45; limrht = 315;
% Elastic Modulus
EMod = 100000;
% Yield force & plastic modulus
% YBar = 1000000; Ydf = 100000; Ydb = 10000000; 
YBar = 5000; Ydf = 0.003; Ydb = 0.03; 
PMod = 100000; pl_mod_fold = 0.0033; pl_mod_bend = 0.033;
% Left and right limits for the linear range of rotational stiffness
[Node,Panel]=ConfigMorph(sec_hor,sec_vert,a,b,alfa1,alfa2,beta,sgngamma2); 

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
indp = 1;
Supp = [];
for hor_id=1:(2*sec_hor+1)
    if hor_id==sec_hor+1
        Supp = [Supp;
                (2*sec_vert+1)*hor_id, 1, 1, 0;];
    elseif mod(hor_id,2)==1
        Supp = [Supp;
                (2*sec_vert+1)*hor_id, 0, 1, 1;];
    else
        Supp = [Supp;
            (2*sec_vert+1)*hor_id, 0, 1, 0;];
    end
end
for vert_id=1:(2*sec_vert+1)
    Supp = [Supp;
            sec_hor*(2*sec_vert+1)+vert_id, 1, 0, 0];
end
% 黄色固定z
for hor_id=1:2:(2*sec_hor+1)
    for vec_id=1:2:(2*sec_vert-1) 
        Supp = [Supp;
                (2*sec_vert+1)*(hor_id-1)+vec_id, 0, 0, 1;];
    end
end
Load = [];
for hor_id=2:2:2*sec_hor
        Load = [Load;
                (2*sec_vert+1)*(hor_id-1)+1, 0, 9, 0;];
end
% Supp(1,2)=1;
% Load = [];
% for hor_id=1:2*sec_hor+1
%     Load = [Load;
%             (2*sec_vert+1)*(hor_id-1)+1, 0, 2, 0;];
% end

% Simulation options using the N4B5 model
DIRC = [0, 1, 0];
% 加载方向
LoadSize = 0.02;
% CycleIcrm 是循环加载中一个加载/一个卸载的步骤数
CycleIcrm = 500;
% TotalIcr 是循环加载的总数
TotalIcr = 500;
AnalyInputOpt = struct(...
    'ModelType','N4B5',...
    'MaterCalib','manual',... 
    'MaxIcr', TotalIcr, ...
    'Abar', Abar,...
    'Kb',Kb,...
    'Kf',Kf,...
    'RotSprBend', @(he,h0,Kb,L0)EnhancedLinear(he,h0,Kb,L0,30,330),...
    'RotSprFold', @(he,h0,Kf,L0)EnhancedLinear(he,h0,Kf,L0,30,330),...
    'LoadType','Displacement',...    % Displacement load
    'DispStep', TotalIcr);

    % 'AdaptiveLoad', @(Node,U,icrm)LoadFunMorph(Node,U,icrm,LoadSize,DIRC,CycleIcrm,sec_hor,sec_vert), ...

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
instdof = [sec_hor*(2*sec_vert+1)+1,2];
interv = 1; endicrm = size(Uhis,2);
% Animation monitoring node-wise change
VIntensityDataInten = zeros(size(truss.Node,1),size(Uhis,2));
IntensityDataM = bsxfun(@times,STAT.bar.Sx,truss.A);
for k = 1:size(Uhis,2)
    IntensityDataIntenk = sparse(truss.Bars(:,1),truss.Bars(:,2),abs(IntensityDataM(:,k)),size(truss.Node,1),size(truss.Node,1));
    VIntensityDataInten(:,k) = sum((IntensityDataIntenk+IntensityDataIntenk'),2); 
end
VisualFold(Uhis(:,1:interv:endicrm),truss,angles,Fhis(1:interv:endicrm,:),instdof,'IntensityMap','Vertex','IntensityData',VIntensityDataInten)
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