clear
%lecture des données 
%[nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient] = lireFichier('instanceExemple.dat');
load('variables.mat')
%% 

%exemple d'appel à la résolution du modèle
[solution, fval] = optimProd(1,nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient);



%%%% PROGRAMMATION DES MODELES (à compléter)%%%%%%%%%%%%%%%
function [solution, fval] = optimProd(modele, nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient)

    %fval=[] ; % ligne à enlever quand un premier modèle sera codé       
    %solution=[] ;% ligne à enlever quand un premier modèle sera codé
    T=30;
    if modele==1
        %C'est pour la premiere partie du sujet
        problem=optimproblem("ObjectiveSense","minimize");
        %Variables de décisions
        x=optimvar('x',nbProduits,T,LowerBound=0);
        s=optimvar('s',nbProduits,T,LowerBound=0);
        y=optimvar('y',nbProduits,nbClients,T);
        %Fonction Objective
        coutStockage=0;
        for i  = 1:nbProduits
            for t = 1:T
                coutStockage = coutStockage+ coutStockUsine(i)*s(i,t);
            end
        end
        coutPenaliteLivraison=0;
        for i = 1:nbProduits
            for j = 1 : nbClients
                for t = 1:T
                    coutPenaliteLivraison=coutPenaliteLivraison+ ((t<a(j))*(a(j)-t)*penalite(j) + (t>b(j))*(t-b(j))*penalite(j))*y(i,j,t);
                end
            end
        end
        funcObjective=coutStockage+coutPenaliteLivraison;
        problem.Objective=funcObjective;
        %Contraintes
        const1=optimconstr(nbClients,T);
        for i = 1:nbProduits
            for t = 1:T
                const1(i,t)= x(i,t)<=capaProd(i);
            end
        end
        problem.Constraints.c1=const1; %Capacité de production
        %problem.Constraints.c1= x<=capaProd; %Capacité de production
        const2=optimconstr(nbClients,T);
        for i = 1 : nbProduits
            for t = 1:T
                if t>1
                    somme=0;
                    for j = 1:nbClients
                        somme=somme+y(i,j,t);
                    end
                end
            end
        end
        problem.Constraints.c2= s == tabCont1;%Equilibre stock
        
        tabCont2=zeros(nbProduits,nbClients);
        for i = 1:nbProduits
            for j=1:nbClients
                tabCont2(i,j)=0;
                for t=1:T
                    tabCont2(i,j)=tabCont2(i,j)+y(i,j,t);
                end
            end
        end
        problem.Constraints.c3=demande == tabCont2;%Satisfaction demande
        tabCont3=zeros(T);
        for t=1:T
            tabCont3(t)=0;
            for i = 1:nbProduits
                for j =1:nbClients
                    tabCont3(t)=tabCont3(t)+y(i,j,t)
                end
            end
        end
        problem.Constraints.c4= tabCont3 <=capaCrossdock;%Capacité entrepôt
        show(problem)
        [sol,fval,exitFlag]=solve(problem,'Solver','linprog');
        if exitFlag>=0
            fprintf('Valeur maximale dans le sac : %d\n', fval);
        end
    elseif modele==2
        %TODO : compléter avec le code de IP1
        fprintf("Pour l'instant, le modèle IP1 n'est pas codé \n"); % ligne à enlever 
    elseif modele==3 
         %TODO : compléter avec le code de IP2
        fprintf("Pour l'instant, le modèle IP2 n'est pas codé \n"); % ligne à enlever 
    else 
        fprintf("Le paramètre modele devrait valoir 1, 2 ou, 3 \n ")
    end
    
end


%%% A compléter
function plotOptim(nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient)
    %TODO : à compléter
     fprintf("Pour l'instant, la fonction n'est pas codée \n"); % ligne à enlever      
end

%%%%%%%FONCTION DE PARSAGE (ne pas modifier)%%%%%%%%
function [nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient]=lireFichier(filename)
% lecture du fichier de données
instanceParameters = fileread(filename);
% suppression des éventuels commentaires
instanceParameters = regexprep(instanceParameters, '/\*.*?\*/', '');
% évaluation des paramètres
eval(instanceParameters);
end


