clear
%lecture des données 
[nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient] = lireFichier('instance1.dat');
%load('variables2.mat')
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
    problem = optimproblem("ObjectiveSense", "minimize");

    % Variables de décision
    x = optimvar('x', nbProduits, T, "LowerBound", 0); % Production
    s = optimvar('s', nbProduits, T, "LowerBound", 0); % Stock
    y = optimvar('y', nbProduits, nbClients, T, "LowerBound", 0); % Livraison

    % Fonction objectif
    coutStockage = sum(coutStockUsine * sum(s, 2));
    coutPenalite = 0;
    for i = 1:nbProduits
        for j = 1:nbClients
            for t = 1:T
                avance = max(0, a(j) - t);
                retard = max(0, t - b(j));
                coutPenalite = coutPenalite + (avance * penalite(j) + retard * penalite(j)) * y(i, j, t);
            end
        end
    end
    problem.Objective = coutStockage + coutPenalite;

    % Contraintes
    % 1. Capacité de production
    for i = 1:nbProduits
        for t = 1:T
            problem.Constraints.("capacite_"+i+"_"+t) = x(i, t) <= capaProd(i);
        end
    end

    % 2. Équilibre des stocks
    for i = 1:nbProduits
        for t = 1:T-1
            problem.Constraints.("equilibre_"+i+"_"+t) = ...
                s(i, t+1) == s(i, t) + x(i, t+1) - sum(y(i, :, t+1));
        end
    end

    % 3. Satisfaction des clients
    for i = 1:nbProduits
        for j = 1:nbClients
            problem.Constraints.("satisfaction_"+i+"_"+j) = ...
                sum(y(i, j, :)) == demande(i, j);
        end
    end

    % 4. Capacité de l'entrepôt central
    for t = 1:T
        problem.Constraints.("capacite_entrepot_"+t) = ...
            sum(sum(y(:, :, t))) <= capaCrossdock;
    end

    % Résolution
    [solution, fval] = solve(problem);
    fprintf("Valeur objective : %f\n", fval);
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


