function [F] = LoadFunKreslingYasuda(Node,U,icrm)
%% Compute current configuration 
Nodenw = Node;
Nodenw(:,1) = Node(:,1)+U(1:3:end);
Nodenw(:,2) = Node(:,2)+U(2:3:end);
Nodenw(:,3) = Node(:,3)+U(3:3:end);
DIRC = [0, 0, 1];
StepLength = 0.00005;
%% Define Load here based on icrm and Nodenw
if icrm<=0 
    error('Wrong increment!'); 
elseif icrm<=50
    Load = [7, StepLength*DIRC/norm(DIRC);
            8, StepLength*DIRC/norm(DIRC);
            9, StepLength*DIRC/norm(DIRC);
            10, StepLength*DIRC/norm(DIRC);
            11, StepLength*DIRC/norm(DIRC);
            12, StepLength*DIRC/norm(DIRC);];
elseif icrm<=450
    Load = [7, -StepLength*DIRC/norm(DIRC);
            8, -StepLength*DIRC/norm(DIRC);
            9, -StepLength*DIRC/norm(DIRC);
            10, -StepLength*DIRC/norm(DIRC);
            11, -StepLength*DIRC/norm(DIRC);
            12, -StepLength*DIRC/norm(DIRC);];
elseif icrm<=850
    Load = [7, StepLength*DIRC/norm(DIRC);
            8, StepLength*DIRC/norm(DIRC);
            9, StepLength*DIRC/norm(DIRC);
            10, StepLength*DIRC/norm(DIRC);
            11, StepLength*DIRC/norm(DIRC);
            12, StepLength*DIRC/norm(DIRC);];
elseif icrm<=1250
    Load = [7, -StepLength*DIRC/norm(DIRC);
            8, -StepLength*DIRC/norm(DIRC);
            9, -StepLength*DIRC/norm(DIRC);
            10, -StepLength*DIRC/norm(DIRC);
            11, -StepLength*DIRC/norm(DIRC);
            12, -StepLength*DIRC/norm(DIRC);];
elseif icrm<=1650
    Load = [7, StepLength*DIRC/norm(DIRC);
            8, StepLength*DIRC/norm(DIRC);
            9, StepLength*DIRC/norm(DIRC);
            10, StepLength*DIRC/norm(DIRC);
            11, StepLength*DIRC/norm(DIRC);
            12, StepLength*DIRC/norm(DIRC);];
elseif icrm<=2050
    Load = [7, -StepLength*DIRC/norm(DIRC);
            8, -StepLength*DIRC/norm(DIRC);
            9, -StepLength*DIRC/norm(DIRC);
            10, -StepLength*DIRC/norm(DIRC);
            11, -StepLength*DIRC/norm(DIRC);
            12, -StepLength*DIRC/norm(DIRC);];
elseif icrm<=2450
    Load = [7, StepLength*DIRC/norm(DIRC);
            8, StepLength*DIRC/norm(DIRC);
            9, StepLength*DIRC/norm(DIRC);
            10, StepLength*DIRC/norm(DIRC);
            11, StepLength*DIRC/norm(DIRC);
            12, StepLength*DIRC/norm(DIRC);];
elseif icrm<=2850
    Load = [7, -StepLength*DIRC/norm(DIRC);
            8, -StepLength*DIRC/norm(DIRC);
            9, -StepLength*DIRC/norm(DIRC);
            10, -StepLength*DIRC/norm(DIRC);
            11, -StepLength*DIRC/norm(DIRC);
            12, -StepLength*DIRC/norm(DIRC);];
elseif icrm<=3250
    Load = [7, StepLength*DIRC/norm(DIRC);
            8, StepLength*DIRC/norm(DIRC);
            9, StepLength*DIRC/norm(DIRC);
            10, StepLength*DIRC/norm(DIRC);
            11, StepLength*DIRC/norm(DIRC);
            12, StepLength*DIRC/norm(DIRC);];
elseif icrm<=3650
    Load = [7, -StepLength*DIRC/norm(DIRC);
            8, -StepLength*DIRC/norm(DIRC);
            9, -StepLength*DIRC/norm(DIRC);
            10, -StepLength*DIRC/norm(DIRC);
            11, -StepLength*DIRC/norm(DIRC);
            12, -StepLength*DIRC/norm(DIRC);];
elseif icrm<=4050
    Load = [7, StepLength*DIRC/norm(DIRC);
            8, StepLength*DIRC/norm(DIRC);
            9, StepLength*DIRC/norm(DIRC);
            10, StepLength*DIRC/norm(DIRC);
            11, StepLength*DIRC/norm(DIRC);
            12, StepLength*DIRC/norm(DIRC);];
elseif icrm<=4450
    Load = [7, -StepLength*DIRC/norm(DIRC);
            8, -StepLength*DIRC/norm(DIRC);
            9, -StepLength*DIRC/norm(DIRC);
            10, -StepLength*DIRC/norm(DIRC);
            11, -StepLength*DIRC/norm(DIRC);
            12, -StepLength*DIRC/norm(DIRC);];
else
    Load = [7, StepLength*DIRC/norm(DIRC);
            8, StepLength*DIRC/norm(DIRC);
            9, StepLength*DIRC/norm(DIRC);
            10, StepLength*DIRC/norm(DIRC);
            11, StepLength*DIRC/norm(DIRC);
            12, StepLength*DIRC/norm(DIRC);];
end

%% Wrap up Load info to force/displacement vector F
m = size(Node,1);
F = zeros(3*m,1);
indp = Load(:,1);
F(3*indp-2) = Load(:,2); 
F(3*indp-1) = Load(:,3); 
F(3*indp) = Load(:,4);
