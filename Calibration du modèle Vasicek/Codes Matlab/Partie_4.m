%%TP3 - Vasicek Model
%%Part 4 - Calibration to historical dates

function[]=Partie_4()
%% Valeurs initiales

eta = 0.6;
sigma = 0.08 ; 
gamma = 4;
T = 5;

N = 1000;
t = linspace(0,T,N); 

r0 = 0;
r(1) = r0;
r_estime(1) = r0;
r = linspace(0,T,N);
r_estime = linspace(0,T,N);


%% Calculer la variance D², gamma, eta, sigma

%Fonction pour le calcul de la variance D²
function[D] = Calcul_Variance_D2(r,a,b)
    D = 0;
    D1 = 0;
    D2 = 0;
    N = length(r)-1;
    for i = 1:N
        D1 = D1+(r(i+1)-(a*r(i)+b))^2;
        D2 = D2+(r(i+1)-(a*r(i)+b));
    end
    D1 = D1/N;
    D2 = D2/N;
    D = D1-D2^2;
    D = sqrt(D);
end


%Fonction pour le calcul de gamma
function[gamma] = Calcul_Gamma(a,dt)
    gamma = - log(a)/dt;
end

%Fonction pour le calcul de eta
function[eta] = Calcul_Eta(gamma,a,b)
    eta = gamma * b/(1-a);
end

%Fonction pour le calcul de sigma
function[sigma] = Calcul_Sigma(D,a,dt)
    sigma = D * sqrt( (- 2 * log(a))/(dt * (1-a^2)));
end


%% Simulation des dates de marché avec les formules théoriques

%----------------------Construction du vecteur r--------------------------
function[r,dt] = Simulation_r(t,gamma,eta,sigma)
    
    r(1) = r0;
    r = linspace(0,T,N);

    for i = 1:N-1
        dt= t(i+1) - t(i);
        u = randn;
        r(i+1) = r(i)*exp(-gamma*dt)+(eta/gamma)*(1-exp(-gamma*dt))+sigma*sqrt((1-exp(-2*gamma*dt))/(2*gamma))*u;
    end
    
end

%Appel de la fonction
[r,dt] = Simulation_r(t,gamma,eta,sigma);



%---Fonction permettant de réaliser le tracé de r en fonction du temps----
function[r] = Graphique_r_Temps(t,r)
    
    %Tracé du graphique de r en fonction du temps
    figure;
    plot(t,r)
    title("Evolution de r en fonction du temps")
    xlabel("Temps")
    ylabel("r")
    
end

%Appel de la fonction
Graphique_r_Temps(t,r)



%-----------Fonction permettant de réaliser le tracé de r(i)-------------
%----------------------------en fonction de r(i+1) ----------------------
function[r] = Graphique_r(r)
    
    %Tracé du graphique r(i) en fonction de r(i+1)
    figure;
    hold on;
    
    for j = 1:N-1
        plot(r(j),r(j+1),"g+");
    end
    
    xlabel("r(i)")
    ylabel("r(i+1)")
    title("Graphique de r(i) en fonction de r(i+1)")
    
end

%Appel de la fonction
Graphique_r(r)



%% Calibration et estimation des paramètres

function[r] = Calibration_Historical_Dates(r,dt)
        
    %Graphique des points r(i) r(i+1)
    Graphique_r(r);
    
    
    %Realisation de la régression linéaire
    x = ones(length(r)-1,2);
    for i = 1: (length(r)-1)
        x(i,1) = r(i);
    end
    
    
    %Détermination des paramètres et tracé de la droite estimée
    y = r(2:length(r))';
    coefficient_regression = x\y;
    a_estime = coefficient_regression(1);
    b_estime = coefficient_regression(2);
    droite_estime = x(:,1) * a_estime + b_estime;
    plot(x(:,1),droite_estime,'blue');
    
    
    %Détermination des paramètres et tracé de la droite théorique
    a_theorique = exp(-gamma*dt);
    b_theorique = (eta/gamma)*(1-exp(-gamma*dt));
    droite_theorique = a_theorique*r+b_theorique;
    plot(r,droite_theorique,'magenta');
    
    
    %Ajout des labels
    xlabel('r')
    ylabel('y=ax+b')
    title('Régression linéaire sur les points de marché selon les coefficients utilisés')
    
    
    %Calcul des paramètres estimés sigma, eta et gamma à partir des
    %formules théoriques mais en prenant les paramètres estimés et trouvés
    %précedemment lors de la régression linéaire pour la droite estimée
    D = Calcul_Variance_D2(r,a_estime,b_estime); 
    gamma_estime = Calcul_Gamma(a_estime,dt);
    eta_estime = Calcul_Eta(gamma_estime, a_estime, b_estime); 
    sigma_estime = Calcul_Sigma(D,a_estime,dt);
    
    
    %Comparaison des paramètres théoriques et estimés
    disp('-----------------------------------------')
    disp(['La valeur de a théorique pour la régression linéaire y=ax+b = ', num2str(a_theorique)])
    disp(['La valeur de a estimée de la régression linéaire y=ax+b = ', num2str(a_estime)])
    disp('-----------------------------------------')
    disp(['La valeur de b théorique pour la régression linéaire y=ax+b = ', num2str(b_theorique)])
    disp(['La valeur de b estimée de la régression linéaire y=ax+b = ', num2str(b_estime)])
    disp('-----------------------------------------')
    disp(['La valeur de la variance D^2 est : ',num2str(D^2)])
    disp('-----------------------------------------')
    disp(['Gamma théorique = ',num2str(gamma)])
    disp(['Gamma estimé = ',num2str(gamma_estime)])
    disp('-----------------------------------------')
    disp(['Eta théorique = ',num2str(eta)])
    disp(['Eta estimé = ',num2str(eta_estime)])
    disp('-----------------------------------------')
    disp(['Sigma théorique = ',num2str(sigma)])
    disp(['Sigma estimé = ',num2str(sigma_estime)])
    disp('-----------------------------------------')
    
end

Calibration_Historical_Dates(r,dt);


end