%calculate Loads 1To1 
%INPUTs:
% nNodes - number of nodes
% Links - matrix containing all Links
% T - flow bandwidth information
% sP - k must available routing paths for each taffic flow
% alternativePaths - k alternative paths for each traffic flow
% sol - indication of each disjount pair we are using
%
% Note: a K must available routing path together with a k alternative
%       path form a disjount pair

function Loads=calculateLoads1To1(nNodes,Links,T,sP,alternativePaths, sol)
    nFlows = size(sP,2);
    nLinks = size(Links,1);
    %link loads using most available path
    Loads= calculateLinkLoads(nNodes,Links,T,sP,sol);
    
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
        sol = zeros(1,nFlows);
        for i = 1:nFlows
            sol(i) = ~isempty(auxSP{i});
        end
        auxLoads= calculateLinkLoads(nNodes,Links,T,auxSP,sol);
        
        Loads1To1(:,3:4) = max(Loads1To1(:,3:4), auxLoads(:,3:4)); 
    end
end