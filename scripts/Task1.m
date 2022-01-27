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

%1.a.
fprintf("\n---------1.a.---------\n");

% Compute up to n paths for each flow:
n= inf;
[sP, nSP]= calculatePaths(L,T,n);

for flow = 1:nFlows
    numberOfPaths = nSP(flow);
    
    %print number of paths for each flow
    fprintf('\nFlow %d number of paths: %d \n ', flow, numberOfPaths);
end

fprintf("\n---------1.b.---------\n");
%1.b.i
%Optimization algorithm resorting to the random strategy:
t= tic;
bestLoad= inf;
worstLoad = -inf;
sol= zeros(1,nFlows);
allValues= [];
while toc(t)<10
    for i= 1:nFlows
        sol(i)= randi(nSP(i));
    end
    Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
    load= max(max(Loads(:,3:4)));
    allValues= [allValues load];
    if load<bestLoad
        bestSol= sol;
        bestLoad= load;
    end
    if load>worstLoad
        worstLoad = load;
    end
end
figure(1)
hold on
plot(sort(allValues));
fprintf('RANDOM STRATEGY:\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));

%1.b.ii
%Optimization algorithm resorting to the random strategy: (10 shortest
%paths)
t= tic;
bestLoad= inf;
worstLoad = -inf;
sol= zeros(1,nFlows);
allValues= [];
while toc(t)<10
    for i= 1:nFlows
        sol(i)= randi(10);
    end
    Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
    load= max(max(Loads(:,3:4)));
    allValues= [allValues load];
    if load<bestLoad
        bestSol= sol;
        bestLoad= load;
    end
    if load>worstLoad
        worstLoad = load;
    end
end
plot(sort(allValues));
fprintf('RANDOM STRATEGY (10 shortest paths):\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));

%1.b.iii
%Optimization algorithm resorting to the random strategy: (5 shortest
%paths)
t= tic;
bestLoad= inf;
worstLoad = -inf;
sol= zeros(1,nFlows);
allValues= [];
while toc(t)<10
    for i= 1:nFlows
        sol(i)= randi(5);
    end
    Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
    load= max(max(Loads(:,3:4)));
    allValues= [allValues load];
    if load<bestLoad
        bestSol= sol;
        bestLoad= load;
    end
    if load>worstLoad
        worstLoad = load;
    end
end
plot(sort(allValues));
legend('Random strategy','Random strategy (10 shortest paths)','Random strategy (5 shortest paths)', 'Location', 'northwest');

fprintf('RANDOM STRATEGY (5 shortest paths):\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));

fprintf("\n---------1.c.---------\n");
%1.c.i
%Optimization algorithm with greedy randomized:
t= tic;
bestLoad= inf;
worstLoad = -inf;
allValues= [];
while toc(t)<10
    ax2= randperm(nFlows);
    sol= zeros(1,nFlows);
    for i= ax2
        k_best= 0;
        best= inf;
        for k= 1:nSP(i)
            sol(i)= k;
            Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
            load= max(max(Loads(:,3:4)));
            if load<best
                k_best= k;
                best= load;
            end
        end
        sol(i)= k_best;
    end
    load= best;
    allValues= [allValues load];
    if load<bestLoad
        bestSol= sol;
        bestLoad= load;
    end
    if load>worstLoad
        worstLoad = load;
    end
end
figure(2)
hold on
plot(sort(allValues));
fprintf('GREEDY RANDOMIZED:\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));

%1.c.ii
%Optimization algorithm with greedy randomized (10 shortest paths):
t= tic;
bestLoad= inf;
worstLoad = -inf;
allValues= [];
while toc(t)<10
    ax2= randperm(nFlows);
    sol= zeros(1,nFlows);
    for i= ax2
        k_best= 0;
        best= inf;
        for k= 1:10 %<--- 1 to 10
            sol(i)= k;
            Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
            load= max(max(Loads(:,3:4)));
            if load<best
                k_best= k;
                best= load;
            end 
        end
        sol(i)= k_best;
    end
    load= best;
    allValues= [allValues load];
    if load<bestLoad
        bestSol= sol;
        bestLoad= load;
    end
    if load>worstLoad
        worstLoad = load;
    end
end

plot(sort(allValues));
fprintf('GREEDY RANDOMIZED (10 shortest paths):\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));

%1.c.iii
%Optimization algorithm with greedy randomized (5 shortest paths):
t= tic;
bestLoad= inf;
worstLoad = -inf;
allValues= [];
while toc(t)<10
    ax2= randperm(nFlows);
    sol= zeros(1,nFlows);
    for i= ax2
        k_best= 0;
        best= inf;
        for k= 1:5    %<--- 1 to 5
            sol(i)= k;
            Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
            load= max(max(Loads(:,3:4)));
            if load<best
                k_best= k;
                best= load;
            end  
        end
        sol(i)= k_best;
    end
    load= best;
    allValues= [allValues load];
    if load<bestLoad
        bestSol= sol;
        bestLoad= load;
    end
    if load>worstLoad
        worstLoad = load;
    end
end
plot(sort(allValues));
legend('Greedy Randomized','Greedy Randomized (10 shortest paths)','Greedy Randomized (5 shortest paths)', 'Location', 'southeast');

fprintf('GREEDY RANDOMIZED (5 shortest paths):\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));

fprintf("\n---------1.d.---------\n");

%1.d.i
%Optimization algorithm with multi start hill climbing:
t= tic;
bestLoad= inf;
worstLoad= -inf;
allValues= [];
while toc(t)<10
    
    %GREEDY RANDOMIZED:
    ax2= randperm(nFlows);
    sol= zeros(1,nFlows);
    for i= ax2
        k_best= 0;
        best= inf;
        for k= 1:nSP(i)
            sol(i)= k;
            Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
            load= max(max(Loads(:,3:4)));
            if load<best
                k_best= k;
                best= load;
            end
        end
        sol(i)= k_best;
    end
    load= best;
    
    %HILL CLIMBING:
    continuar= true;
    while continuar
        i_best= 0; %fluxo
        k_best= 0; %percurso do fluxo
        best= load;
        for i= 1:nFlows
            for k= 1:nSP(i)
                if k~=sol(i)
                    aux= sol(i);
                    sol(i)= k;
                    Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
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
    if load>worstLoad
        worstLoad = load;
    end
end
figure(3)
hold on
plot(sort(allValues));

fprintf('MULTI START HILL CLIMBING:\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));


%1.d.ii
%Optimization algorithm with multi start hill climbing (10 shortest paths):
t= tic;
bestLoad= inf;
worstLoad= -inf;
allValues= [];
while toc(t)<10
    
    %GREEDY RANDOMIZED:
    ax2= randperm(nFlows);
    sol= zeros(1,nFlows);
    for i= ax2
        k_best= 0;
        best= inf;
        for k= 1:10   % <--- 1 to 10
            sol(i)= k;
            Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
            load= max(max(Loads(:,3:4)));
            if load<best
                k_best= k;
                best= load;
            end
        end
        sol(i)= k_best;
    end
    load= best;
    
    %HILL CLIMBING:
    continuar= true;
    while continuar
        i_best= 0; %fluxo
        k_best= 0; %percurso do fluxo
        best= load;
        for i= 1:nFlows
            for k= 1:10 % <--- 1 to 10
                if k~=sol(i)
                    aux= sol(i);
                    sol(i)= k;
                    Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
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
    if load>worstLoad
        worstLoad = load;
    end
end

plot(sort(allValues));

fprintf('MULTI START HILL CLIMBING (10 shortest paths):\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));


%1.d.iii
%Optimization algorithm with multi start hill climbing (5 shortest paths):
t= tic;
bestLoad= inf;
worstLoad= -inf;
allValues= [];
while toc(t)<10
    
    %GREEDY RANDOMIZED:
    ax2= randperm(nFlows);
    sol= zeros(1,nFlows);
    for i= ax2
        k_best= 0;
        best= inf;
        for k= 1:5    % <---- 1 to 5
            sol(i)= k;
            Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
            load= max(max(Loads(:,3:4)));
            if load<best
                k_best= k;
                best= load;
            end
        end
        sol(i)= k_best;
    end
    load= best;
    
    %HILL CLIMBING:
    continuar= true;
    while continuar
        i_best= 0; %fluxo
        k_best= 0; %percurso do fluxo
        best= load;
        for i= 1:nFlows
            for k= 1:5   % <--- 1 to 5
                if k~=sol(i)
                    aux= sol(i);
                    sol(i)= k;
                    Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
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
    if load>worstLoad
        worstLoad = load;
    end
end

plot(sort(allValues));
legend('Hill Climbing', 'Hill Climbing (10 shortest paths)','Hill Climbing (5 shortest paths)', 'Location', 'northwest');

fprintf('MULTI START HILL CLIMBING (5 shortest paths):\n');
fprintf('   Best load = %.2f Gbps\n',bestLoad);
fprintf('   Worst load = %.2f Gbps\n',worstLoad);
fprintf('   No. of solutions = %d\n',length(allValues));
fprintf('   Av. quality of solutions = %.2f Gbps\n',mean(allValues));