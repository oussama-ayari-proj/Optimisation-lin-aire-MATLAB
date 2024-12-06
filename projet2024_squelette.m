clear
%lecture des données 
[nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient] = lireFichier('instanceExemple.dat');
%load('variables2.mat')
%%
% Appel pour la premiére partie
%[solution, fval] = optimProd(1,nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient);
% Appel pour la question 1
[~, ~] = optimProd(1,nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient);
% Appel pour la question 5
[~, ~] = optimProd(2,nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient);
% Appel pour la question 6
[~, ~] = optimProd(3,nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient);
% Appel pour la question 3
%plotOptim(nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient)

% Début fonction Question 2
function T = calculerHorizon(I, F, d, b, M)
    % Étendre jusqu'à la fenêtre maximale de livraison de livraison
    TmaxLivraison = max(b);
    % Calculer le temps minimal requis pour la production
    Tprod = 0;
    for i = 1:I
        Tprod = max(Tprod, sum(d(i, :)) / F(i));
    end
    % Calculer le temps minimal pour gérer toutes les livraisons
    Tentrepot = ceil(sum(sum(d)) / M);
    % Prendre le maximum des trois
    T = max([TmaxLivraison, Tprod, Tentrepot]);
end
% Fin fonction Question 2

%%%% PROGRAMMATION DES MODELES (à compléter)%%%%%%%%%%%%%%%
function [solution, fval] = optimProd(modele, nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient)
    % Appel pour la Question 2 (Calcul de l'horizon optimal T)
    %T=30; % Valeur par défaut de T
    T_optim=calculerHorizon(nbProduits,capaProd,demande,b,capaCrossdock);
    T=T_optim;
    % Fin Appel de la Question 2
    %C'est pour la premiere partie du sujet
    % Début Question 1 (Modélisation du probléme)
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
                avance=0;
                if t<a(j)
                    avance=a(j)-t;
                end
                retard=0;
                if t>b(j)
                    retard=t-b(j);
                end
                coutPenalite = coutPenalite + (avance * penalite(j) + retard * penalite(j)) * y(i, j, t);
            end
        end
    end
    problem.Objective = coutStockage + coutPenalite;

    % Contraintes
    % Capacité de production
    for i = 1:nbProduits
        for t = 1:T
            problem.Constraints.("capacite_"+i+"_"+t) = x(i, t) <= capaProd(i);
        end
    end

    % Équilibre des stocks
    for i = 1:nbProduits
        for t = 1:T-1
            problem.Constraints.("equilibre_"+i+"_"+t) = ...
                s(i, t+1) == s(i, t) + x(i, t+1) - sum(y(i, :, t+1));
        end
    end

    % Satisfaction des clients (On suppose qu'on va satisfaire tous les
    % clients)
    for i = 1:nbProduits
        for j = 1:nbClients
            problem.Constraints.("satisfaction_"+i+"_"+j) = ...
                sum(y(i, j, :)) == demande(i, j);
        end
    end

    % Capacité de l'entrepôt central
    for t = 1:T
        problem.Constraints.("capacite_entrepot_"+t) = ...
            sum(sum(y(:, :, t))) <= capaCrossdock;
    end
    if modele==1
    % Résolution
    [solution, fval] = solve(problem,"Solver","linprog");
    fprintf("Valeur objective du model %d : %f \n",modele,fval);
    % Fin Question 1
    elseif modele==2
        % Début de la question 5
        z = optimvar('z', nbProduits, T, "Type", "integer", "LowerBound", 0, "UpperBound", 1); 
        w = optimvar('w', nbClients, T, "Type", "integer", "LowerBound", 0, "UpperBound", 1);
        % Fonction objectif
        % Cout de transport
        coutTransport = 0;
        for i = 1:nbProduits
            for t = 1:T
                coutTransport = coutTransport + coutCamionUsine(i) * z(i, t);
            end
        end
        for j = 1:nbClients
            for t = 1:T
                coutTransport = coutTransport + coutCamionClient(j) * w(j, t);
            end
        end

        problem.Objective = coutStockage + coutPenalite + coutTransport;
        % Contraintes
        % Disponibilité des Camions (usine entrepôt)
        big_num=1e6;
        for i = 1:nbProduits
            for t = 1:T
                problem.Constraints.("dispo_camion_usine_"+i+"_"+t) = ...
                sum(y(i, :, t))<=z(i,t)*big_num;
                %z(i,t) >= sum(y(i, :, t)) / big_num;
            end
        end
        % Disponilibilité des Camions (entrepôt client)
        for j = 1:nbClients
            for t = 1:T
                problem.Constraints.("dispo_camion_client_"+j+"_"+t) = ...
                sum(y(:, j, t)) <= w(j,t)*big_num;
                %w(j,t) >= sum(y(:, j, t)) / big_num;
            end
        end
        % Résolution
        options = optimoptions('intlinprog','Display','none','ConstraintTolerance',1e-8,'AbsoluteGapTolerance',1e-8);
        [solution, fval] = solve(problem,Options=options);
        fprintf("Valeur objective du model %d : %f \n",modele,fval);
        
        % Fin de la question 5
    elseif modele==3 
        % Début de la question 6
        z = optimvar('z', nbProduits,nbClients, T, "Type", "integer", "LowerBound", 0, "UpperBound", 1); 
        w = optimvar('w', nbClients,nbProduits, T, "Type", "integer", "LowerBound", 0, "UpperBound", 1);
        % Fonction objectif
        % Cout de transport
        coutTransport = 0;
        for i = 1:nbProduits
            for j=1:nbClients
                for t = 1:T
                    coutTransport = coutTransport + coutCamionUsine(i) * z(i,j, t);
                end
            end
        end
        for j = 1:nbClients
            for i=1:nbProduits
                for t = 1:T
                    coutTransport = coutTransport + coutCamionClient(j) * w(j,i, t);
                end
            end
        end

        problem.Objective = coutStockage + coutPenalite + coutTransport;
        % Contraintes
        % Disponibilité des Camions (usine entrepôt)
        big_num=1e8;
        for i = 1:nbProduits
            for j = 1:nbClients
                for t = 1:T
                    problem.Constraints.("dispo_camion_usine_"+i+"_"+j+"_"+t) = ...
                    y(i, j, t)<= z(i,j, t)*big_num;    
                    %z(i, t) >= y(i, j, t) / big_num;
                end
            end
        end
        % Disponilibilité des Camions (entrepôt client)
        for i = 1:nbProduits
            for j = 1:nbClients
                for t = 1:T
                    problem.Constraints.("dispo_camion_client_"+i+"_"+j+"_"+t) = ...
                    y(i, j, t) <= w(j,i, t)*big_num;   
                    %w(j, t) >= y(i, j, t) / big_num;
                end
            end
        end
        % Résolution
        options = optimoptions('intlinprog','Display','none','ConstraintTolerance',1e-8,'AbsoluteGapTolerance',1e-8);
        [solution, fval] = solve(problem,Options=options);
        fprintf("Valeur objective du model %d : %f \n",modele,fval);
        
        % Fin de la question 6
    else 
        fprintf("Le paramètre modele devrait valoir 1, 2 ou, 3 \n ")
    end
    
end


%%% A compléter
% Début Question 3
function plotOptim(nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient)
    
    valeurs_M=100:20:300;
    valeurs_objectives=zeros(size(valeurs_M));
    for i= 1:length(valeurs_M)
        M=valeurs_M(i);
        [~,fval]=optimProd(1,nbProduits, nbClients, capaProd, M, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient);
        valeurs_objectives(i)=fval;
    end
    figure;
    plot(valeurs_M, valeurs_objectives, '-o');
    xlabel('Capacité de l''entrepôt (M)');
    ylabel('Valeur de la fonction objectif');
    title('Impact de la capacité de l''entrepôt sur la fonction objectif');
    grid on;
end
% Fin Question 3

%%%%%%%FONCTION DE PARSAGE (ne pas modifier)%%%%%%%%
function [nbProduits, nbClients, capaProd, capaCrossdock, demande, a, b, penalite, coutStockUsine, coutCamionUsine, coutCamionClient]=lireFichier(filename)
% lecture du fichier de données
instanceParameters = fileread(filename);
% suppression des éventuels commentaires
instanceParameters = regexprep(instanceParameters, '/\*.*?\*/', '');
% évaluation des paramètres
eval(instanceParameters);
end