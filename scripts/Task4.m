clear all;
close all;

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

% 4.a. 
%   For each flow, compute 10 pairs of link disjoint paths in the following way. 
%   With a kshortest path algorithm, first compute the k = 10 most available routing paths provided by
%   the network to each traffic flow. Then, compute the most available path which is link
%   disjoint with each of the k previous paths.
fprintf("\n---------4.a.---------\n");

%calculate 10 most available path of each flow
n = 10;
[sP, nSP]= calculatePaths(Alog,T,n);

%compute the most available path disjoint with each of k must available
%paths
disjointPaths = cell(1,nFlows);

for flow = 1:nFlows
    for k = 1:length(sP{flow})
        path = sP{flow}{k};

        tmpAlog = Alog;
        %remove all links of k most available path 
        for node = 2:size(sP{flow}{1},2)
            tmpAlog(sP{flow}{k}(node), sP{flow}{k}(node-1)) = inf;
            tmpAlog(sP{flow}{k}(node-1), sP{flow}{k}(node)) = inf;
        end
        
        %recalculate paths for each flow
        [tmpSP, tmpNSP] = calculatePaths(tmpAlog,T,1);
    
         if ~isempty(tmpSP{flow})
            disjointPaths{flow}{k} = tmpSP{flow}{1};
         end
    end
end

%print all pairs
for flow = 1:nFlows
    fprintf('\n-----Flow %d-----\n\n', flow);
    for k = 1:length(sP{flow})
        fprintf("%d most available path: ", k);
        fprintf('[ ');
        fprintf("%g ",sP{flow}{k});
        fprintf(']\n');
        fprintf("Disjoint path: ");
        fprintf('[ ');
        fprintf("%g ",disjointPaths{flow}{k});
        fprintf(']\n\n');
    end
end

% 4.b.
%   Develop a multi start hill climbing algorithm for this optimization problem using the 10
%   pairs of link disjoint paths computed in 4.a for each flow. Run the algorithm during 30
%   seconds. Present the pair of routing paths of each flow (and its availability) and the
%   average service availability of the best solution. Present the highest required bandwidth
%   value among all links. Compare this solution with the one in 3.d and take all possible
%   conclusions.
fprintf("\n---------4.b.---------\n");

time = 30;

%-----Greedy randomized algorithm-----
t= tic;
bestLoad= inf;
sol= zeros(1,nFlows);
allValues= [];
while toc(t)<time
    sol= zeros(1,nFlows);
    for i = randperm(nFlows)
        k_best = 0;
        best = inf; 
        for k = 1:nSP(i) %pair each pair calculate Loads1To1
            sol(i) = k;
            Loads= calculateLoads1To1(nNodes,Links,T,sP,disjointPaths,sol);
            load= max(max(Loads(:,3:4)));
            if load < best
               k_best = k;
               best = load;
            end
        end
        sol(i) = k_best;
    end
    
    allValues= [allValues best];
    if load<bestLoad
        bestSol= sol;
        bestLoad= best;
    end
end
plot(sort(allValues));
hold on

fprintf('\n------Greedy randomized------\n');
printResults(bestLoad, allValues, sP, disjointPaths, sol, A);

%-----Multi start hill climbing algorithm-----
t= tic;
bestLoad= inf;
sol= zeros(1,nFlows);
allValues= [];
while toc(t)<time
    %Optimization algorithm resorting to the greedy randomized strategy:
    sol= zeros(1,nFlows);
    for i = randperm(nFlows)
        k_best = 0;
        best = inf; 
        for k = 1:nSP(i) %pair each pair calculate Loads1To1
            sol(i) = k;
            Loads= calculateLoads1To1(nNodes,Links,T,sP,disjointPaths,sol);
            load= max(max(Loads(:,3:4)));
            if load < best
               k_best = k;
               best = load;
            end
        end
        sol(i) = k_best;
    end
    load = best;
    
    %HILL CLIMBING:
    continuar= true;
    while continuar
        i_best= 0;      %fluxo
        k_best= 0;      %percurso do fluxo
        best= load;
        for i= 1:nFlows
            for k= 1:nSP(i)
                if k~=sol(i)
                    aux= sol(i);
                    sol(i)= k;
                    Loads= calculateLoads1To1(nNodes,Links,T,sP,disjointPaths, sol);
                    load1= max(max(Loads(:,3:4)));
                    if load1<best
                        i_best= i;
                        k_best= k;
                        best= load1;
                    end
                    sol(i)= aux;
                end
            end
        end
        if i_best>0
            sol(i_best)= k_best;
            load= best;
        else
            continuar= false;
        end
    end
    allValues= [allValues load];
    if load<bestLoad
        bestSol= sol;
        bestLoad= load;
    end
end
plot(sort(allValues));
hold on

fprintf('\n------Multi Start Hill Climbing------\n');
printResults(bestLoad, allValues, sP, disjointPaths, sol, A);


%-----Multi start hill climbing algorithm first neighbor-----

t= tic;
bestLoad= inf;
sol= zeros(1,nFlows);
allValues= [];
while toc(t)<time
    %Optimization algorithm resorting to the greedy randomized strategy:
    sol= zeros(1,nFlows);
    for i = randperm(nFlows)
        k_best = 0;
        best = inf; 
        for k = 1:nSP(i) %pair each pair calculate Loads1To1
            sol(i) = k;
            Loads= calculateLoads1To1(nNodes,Links,T,sP,disjointPaths,sol);
            load= max(max(Loads(:,3:4)));
            if load < best
               k_best = k;
               best = load;
            end
        end
        sol(i) = k_best;
    end
    load = best;
    
    %HILL CLIMBING:
    improved = true;
    while improved
        i_best= 0;      %fluxo
        k_best= 0;      %percurso do fluxo
        best= load;

        betterNeighbourFlag = false;
        flow = 1;
        while ~betterNeighbourFlag && flow < nFlows
            nPaths = nSP(flow);
            path = 1;
            while ~betterNeighbourFlag && path < nPaths
                if k~=sol(i)
                    aux= sol(i);
                    sol(i)= k;
                    Loads= calculateLoads1To1(nNodes,Links,T,sP,disjointPaths, sol);
                    load1= max(max(Loads(:,3:4)));
                    if load1<best
                        i_best= i;
                        k_best= k;
                        best= load1;
                        betterNeighbourFlag = true;
                    end
                    sol(i)= aux;
                end
                path = path + 1;
            end
            flow = flow + 1;
        end
        if i_best>0
            sol(i_best)= k_best;
            load= best;
        else
            improved= false;
        end
    end
    allValues= [allValues load];
    if load<bestLoad
        bestSol= sol;
        bestLoad= load;
    end
end
plot(sort(allValues));
legend("Greedy Randomized heuristic","Multi start hill climbing heuristic", ...
       "Multi start hill climbing heuristic first neighbor", "location", "best");

fprintf('\n------Multi Start Hill Climbing First Neighbor------\n');
printResults(bestLoad, allValues, sP, disjointPaths, sol, A);


%display results for 4.b. resolution
function printResults(bestLoad, allValues, sP, disjointPaths, sol, A)
    nFlows = size(sP,2);
    
    pathA1 = ones(1,nFlows);      %must available path availability
    pathA2 = ones(1,nFlows);      %disjoint path availability

    for flow = 1:nFlows
        % Compute flow availability
        for node = 2:size(sP{flow}{1},2)
            availability = A(sP{flow}{1}(node), sP{flow}{1}(node-1)); 
            pathA1(flow) = pathA1(flow) * availability;
        end
    end

    for flow = 1:nFlows
        % Compute flow availability
        for node = 2:size(disjointPaths{flow}{1},2)
            availability = A(disjointPaths{flow}{1}(node), disjointPaths{flow}{1}(node-1)); 
            pathA2(flow) = pathA2(flow) * availability;
        end
    end

    fprintf('   Highest required bandwidth= %.2f Gbps\n',bestLoad);
    fprintf('   No. of solutions = %d\n',length(allValues));
    fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));
    
    fprintf("   Solution = \n");
    for flow = 1:nFlows
        k = sol(flow);
        fprintf("\n     Flow %d pair:\n", flow)
        fprintf('       [ ');
        fprintf("%g ",sP{flow}{k});
        fprintf(']\n');
        fprintf("       Must available path availability: %.5f%%\n", pathA1(flow)*100);
        fprintf('       [ ');
        fprintf("%g ",disjointPaths{flow}{k});
        fprintf(']\n');
        fprintf("       Disjoint path availability: %.5f%%\n", pathA2(flow)*100);
    end
    % ??? 
    meanA = mean((pathA1 + pathA2) ./ 2);
    fprintf("\n???Average availability value among all flows of the service: %.5f%%???\n", meanA*100);                      
end