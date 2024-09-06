function [F] = LoadFunMorph(Node,U,icrm,LoadSize,DIRC,CycleIcrm,sec_hor,sec_vert)
%% Compute current configuration 
Nodenw = Node;
Nodenw(:,1) = Node(:,1)+U(1:3:end);
Nodenw(:,2) = Node(:,2)+U(2:3:end);
Nodenw(:,3) = Node(:,3)+U(3:3:end);
%% Define Load here based on icrm and Nodenw

Cycle_id = fix(icrm/CycleIcrm);
Load = [];
for hor_id=1:2*sec_hor+1
    Load = [Load;
            (2*sec_vert+1)*(hor_id-1)+1, (-1)^Cycle_id*LoadSize*DIRC/norm(DIRC);];
end

%% Wrap up Load info to force/displacement vector F
m = size(Node,1);
F = zeros(3*m,1);
indp = Load(:,1);
F(3*indp-2) = Load(:,2); 
F(3*indp-1) = Load(:,3); 
F(3*indp) = Load(:,4);
