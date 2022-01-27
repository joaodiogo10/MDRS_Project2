%calculate Loads 1To1 
%INPUTs:
% nNodes - number of nodes
% Links - matrix containing all Links
% T - flow bandwidth information
% sP - k must available routing paths for each taffic flow
% disjointPaths - k disjoint paths for each traffic flow
% sol - indication of each disjount pair we are using
%
% Note: a K must available routing path together with a k disjoint
%       path form a disjount pair

function Loads1To1=calculateLoads1To1(nNodes,Links,T,sP,disjointPaths, sol)
    nFlows = size(sP,2);
    nLinks = size(Links,1);
    %link loads using 1º path
    Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
    
    Loads1To1 = Loads;
    %calculate loads in case of single link failure 
    for link = 1:nLinks
        node1= Links(link,1);
        node2= Links(link,2);

        auxSP = cell(1,nFlows);

        for flow = 1:nFlows
            k = sol(flow);
            if(k > 0)
                %check if link is in 1º solution path
                path= sP{flow}{k};
                pathdif= find(path==node1 | path==node2);
                if length(pathdif)<2 || pathdif(2)-pathdif(1)>1
                    %link is not in 1º solution path
                    auxSP{flow}{k} = path;
                elseif ~isempty(disjointPaths{flow}{k})
                    %link is in 1º solution path and we have a disjoint path
                    auxSP{flow}{k} = disjointPaths{flow}{k};
                end
            end
        end
    
        %Calculate link load for a particular link failure
        sol = zeros(1,nFlows);
        for i = 1:nFlows
            sol(i) = ~isempty(auxSP{i});
        end
        auxLoads= calculateLinkLoads(nNodes,Links,T,auxSP,sol);
        
        Loads1To1(:,3:4) = max(Loads1To1(:,3:4), auxLoads(:,3:4)); 
    end
end