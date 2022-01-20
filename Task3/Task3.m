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


nNodes= 10;
nNodes= 9;
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
    fprintf('\nFlow %d:\n', flow);
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
fprintf("\nAverage availability value among all flows of the service: %.5f%%\n", meanA*100);
