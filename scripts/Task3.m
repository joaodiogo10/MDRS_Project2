Nodes= [30 70
       350 40
       550 180
       310 130
       100 170
       540 290
       120 240
       400 310
       220 370
       550 380];
   
Links= [1 2
        1 5
        2 3
        2 4
        3 4
        3 6
        3 8
        4 5
        4 8
        5 7
        6 8
        6 10
        7 8
        7 9
        8 9
        9 10];

T= [1  3  1.0 1.0
    1  4  0.7 0.5
    2  7  2.4 1.5
    3  4  2.4 2.1
    4  9  1.0 2.2
    5  6  1.2 1.5
    5  8  2.1 2.5
    5  9  1.6 1.9
    6 10  1.4 1.6];


nNodes= size(Nodes,1);
nLinks= size(Links,1);
nFlows= size(T,1);

co= Nodes(:,1)+1i*Nodes(:,2);

L= inf(nNodes);    %Square matrix with arc lengths (in Km)
for i=1:nNodes
    L(i,i)= 0;
end

C= zeros(nNodes);  %Square matrix with arc capacities (in Gbps)
for i=1:nLinks
    C(Links(i,1),Links(i,2))= 10;  %Gbps
    C(Links(i,2),Links(i,1))= 10;  %Gbps
    d= abs(co(Links(i,1))-co(Links(i,2)));
    L(Links(i,1),Links(i,2))= d+5; %Km
    L(Links(i,2),Links(i,1))= d+5; %Km 
end
L= round(L);  %Km

%MTBF calculation
MTBF= (450*365*24)./L;
A= MTBF./(MTBF + 24);
A(isnan(A))= 0;

%log transformation
Alog = -log(A);

%3.a. 
%   For each flow, compute one of its routing paths given by the most available path.
fprintf("\n---------3.a.---------\n");

%calculate most available path of each flow
[sP, nSP]= calculatePaths(Alog,T,1);

for flow = 1:nFlows
    path = sP{flow}{1};
    
    %print flow path
    fprintf('\nFlow %d most available path:\n', flow);
    fprintf('[ ');
    fprintf("%g ",path);
    fprintf(']\n');
end


%3.b.   
%   For each flow, compute another routing path given by the most available path which is
%   link disjoint with the previously computed routing path. Compute the availability
%   provided by each pair of routing paths. Present all pairs of routing paths of each flow and
%   their availability. Present also the average service availability (i.e., the average availability
%   value among all flows of the service).
fprintf("\n---------3.b.---------\n");

%calculate most available path of each flow
[sP, nSP]= calculatePaths(Alog,T,1);

alternativePaths = cell(1,nFlows);

pathA1 = ones(1,nFlows);      %must available path availability
pathA2 = ones(1,nFlows);      %alternative path availability

for flow = 1:nFlows   
    %compute availability must available path
    for node = 2:size(sP{flow}{1},2)
        availability = A(sP{flow}{1}(node), sP{flow}{1}(node-1)); 
        pathA1(flow) = pathA1(flow) * availability;
    end

    tmpAlog = Alog;
    %remove all links of most available path 
    for node = 2:size(sP{flow}{1},2)
        tmpAlog(sP{flow}{1}(node), sP{flow}{1}(node-1)) = inf;
        tmpAlog(sP{flow}{1}(node-1), sP{flow}{1}(node)) = inf;
    end
    
    %recalculate path for each flow
    [tmpSP, tmpNSP] = calculatePaths(tmpAlog,T,1);
    
    if ~isempty(tmpSP{flow})
        alternativePaths{flow} = {tmpSP{flow}{1}};

        % Compute flow availability
        for node = 2:size(alternativePaths{flow}{1},2)
            availability = A(alternativePaths{flow}{1}(node), alternativePaths{flow}{1}(node-1)); 
            pathA2(flow) = pathA2(flow) * availability;
        end
    end

    %print must available path and it's availability
    fprintf('\nFlow %d:\n', flow);
    fprintf('   Must available path: ');
    fprintf('[ ');
    fprintf("%g ",sP{flow}{1});
    fprintf(']\n');
    fprintf("   Must available path availability: %.5f%%\n", pathA1(flow)*100);

    %print alternative path and it's availability
    fprintf('   Alternative path: ');
    if isempty(tmpSP{flow})
        fprintf("No alternative disjount path\n\n")
    else
        fprintf('[ ');
        fprintf("%g ",alternativePaths{flow}{1});
        fprintf(']\n');
        fprintf("   Alternative path availability: %.5f%%\n", pathA2(flow)*100);
    end
end

% ??? 
meanA = mean((pathA1 + pathA2) ./ 2);
fprintf("\n???Average availability value among all flows of the service: %.5f%%???\n", meanA*100);


% 3.c.
%   Recall that the capacity of all links is 10 Gbps in each direction. Compute how much
%   bandwidth is required on each direction of each link to support all flows with 1+1
%   protection using the previous computed pairs of link disjoint paths. Compute also the total
%   bandwidth required on all links. Register which links do not have enough capacity. 
fprintf("\n---------3.c.---------\n");

%link loads using most available path
[sP, nSP]= calculatePaths(Alog,T,1);
Loads= calculateLinkLoads(nNodes,Links,T,sP,ones(1,nFlows));

%calculate link loads for alternative paths
%only use flows with an alternative path
sol = zeros(1,nFlows);
for i = 1:nFlows
    sol(i) = ~isempty(alternativePaths{i});
end
newLoads= calculateLinkLoads(nNodes,Links,T,alternativePaths,sol);

Loads1Plus1 = Loads;
Loads1Plus1(:,3:4) = Loads(:,3:4) + newLoads(:,3:4);
totalLoad = sum(sum(Loads1Plus1(:,3:4)));

fprintf("Loads 1 plus 1:\n");
disp(Loads1Plus1);
fprintf("\nTotal bandwidth required on all links: %0.2f\n\n", totalLoad);

fprintf("Links that don't have enough capacity: \n");
for link = 1:nLinks
    if(Loads1Plus1(link,3) > 10 || Loads1Plus1(link,4) > 10)
        fprintf("%d - %d\n", Loads1Plus1(link,1), Loads1Plus1(link,2))
    end
end

% 3.d.
%   Compute how much bandwidth is required on each link to support all flows with 1:1
%   protection using the previous computed pairs of link disjoint paths. Compute also the total
%   bandwidth required on all links. Register which links do not have enough capacity and
%   the highest bandwidth value required among all links.
fprintf("\n---------3.d.---------\n");

%link loads using most available path
[sP, nSP]= calculatePaths(Alog,T,1);
Loads= calculateLinkLoads(nNodes,Links,T,sP,ones(1,nFlows));

Loads1To1 = Loads;
%calculate loads in case of single link failure 
for link = 1:nLinks
    node1= Links(link,1);
    node2= Links(link,2);
    
    auxSP = cell(1,nFlows);
    for flow = 1:nFlows
        %check if link is in must available path
        path= sP{flow}{1};
        pathdif= find(path==node1 | path==node2);
        if length(pathdif)<2 || pathdif(2)-pathdif(1)>1
            %link is not in must available path
            auxSP{flow}{1} = path;
        elseif ~isempty(alternativePaths{flow})
            %link is in must available path and we have an alternative
            auxSP{flow}{1} = alternativePaths{flow}{1};
        end 
    end

    %Calculate link load for a particular link failure
    %only use flows with an alternative path
    sol = zeros(1,nFlows);
    for i = 1:nFlows
        sol(i) = ~isempty(auxSP{i});
    end
    auxLoads= calculateLinkLoads(nNodes,Links,T,auxSP,sol);
    
    Loads1To1(:,3:4) = max(Loads1To1(:,3:4), auxLoads(:,3:4)); 
end

totalLoad = sum(sum(Loads1To1(:,3:4)));
fprintf("Loads 1 to 1:\n");
disp(Loads1To1);
fprintf("\nTotal bandwidth required on all links: %0.2f\n\n", totalLoad);

fprintf("Links that don't have enough capacity: \n");
for link = 1:nLinks
    if(Loads1To1(link,3) > 10 || Loads1To1(link,4) > 10)
        fprintf("%d - %d\n", Loads1Plus1(link,1), Loads1Plus1(link,2))
    end
end

