function[]=Partie_1()
%%PARTIE ACTIF PORTEFEUILLE PANDL VAR 


%% Valeurs des données fixes ou initiales

S(1)=1;
r=0.05;
Sigma=0.5;
T=5;
N=100;
K=1.5;
t=linspace(0,T,N+1);
deltat=T/N; 
Nmc=1000;


%% Simulation du prix de l'actif

function[S]=Prix_Actif()
    S(1)=1;
    for i=1:N
        S(i+1)=S(i)*exp((r-(Sigma^2)/2)*deltat+Sigma*sqrt(deltat)*randn);
    end
end


%Simulation de Nmc chemins correspondants au prix de l’actif en partant toujours de S0

function[]=Simulation_Prix_Actif()
    for k=1:Nmc
        S=Prix_Actif();
        plot(t,S,'b');
        hold on;
        title('Simulation de Nmc chemins correspondants au prix de l actif avec pour départ S0 = 1');
        xlabel('Temps t');
        ylabel('Prix de l actif S');
    end
end

Simulation_Prix_Actif();


%Calcul de l'espérance et la variance de ST

function[Esp,Var] = Estimation_Esperance_Variance_ST()
    CalcEsp=0;
    CalcVar=0;
    for k=1:Nmc
        S=Prix_Actif();
        CalcEsp=CalcEsp+S(N+1);
        CalcVar=CalcVar+S(N+1)^2; 
    end
    Esp=CalcEsp/Nmc;
    Var=(CalcVar/Nmc)-Esp^2;
end

[Esp,Var]=Estimation_Esperance_Variance_ST();
disp('Estimation de');
disp('Esperance de ST = ');
disp(Esp);
disp('Variance de ST = ');
disp(Var);

EspTheorique = exp(r*T) ;
VarTheorique = exp(2*r*T)*(exp(Sigma^2)-1) ;
disp('Valeurs théoriques de');
disp('Esperance de ST = ');
disp(EspTheorique);
disp('Variance de ST = ');
disp(VarTheorique);


%% Simulation du portefeuille de couverture

function[f]=Nx(x)
    f=0.5*(1+erf(x/sqrt(2)));
end

function[f]=d1(t,S)
    f=(log(S/K)+(r+(Sigma^2)/2)*(T-t))/(Sigma*sqrt(T-t));
end

function[f]=d2(t,S)
    f=(log(S/K)+(r-(Sigma^2)/2)*(T-t))/(Sigma*sqrt(T-t));
end

function[f]=ValTheoriqueBS(t,S)
    f=S*Nx(d1(t,S))-K*exp(-r*(T-t))*Nx(d2(t,S));
end

%Calcul de la quantité initiale des actions (ou le nombre de parts de
%sous_jacent) :
disp('A0 =') 
A(1)=Nx(d1(t(1),S(1)));
disp(A(1))

%Valeur du Cash
B(1)=1; %Pour la simulation du portefeuille de couverture

%Valeur du portefeuille de départ : 
P(1)=A(1)*S(1)+B(1);
disp('Portefeuille de départ: P0 = A0S0+B0 =')
disp(P(1))

%Calcul du portefeuille actualisé de départ :
V(1)=ValTheoriqueBS(t(1),S(1));
Pactualise(1)=V(1); 
disp('Portefeuille actualisé de départ: P0actualise =  V0 =')
disp(Pactualise(1))


%Une simulation du portefeuille de couverture et prix de l'option, ratio,
%etc...

function[Pactualise,V,A]=Simulation_Prtf_Couv()
    
    A(1)=Nx(d1(t(1),S(1)));
    P(1)=A(1)*S(1)+B(1);
    V(1)=ValTheoriqueBS(t(1),S(1));
    Pactualise(1)=V(1);
    
    for i=1:N
        S(i+1)=S(i)*exp((r-(Sigma^2)/2)*deltat+Sigma*sqrt(deltat)*randn); %Calcul du prix de l'action
        P(i+1)=A(i)*S(i+1)+(1+r*deltat)*B(i); %Valeur du portefeuille
        A(i+1)=Nx(d1(t(i+1),S(i+1))); %Quantité d'actions qu'il faut avoir pour couvrir une option
        B(i+1)=P(i+1)-A(i+1)*S(i+1); %Valeur du cash
        Pactualise(i+1)=P(i+1)+(V(1)-P(1))*exp(r*t(i+1)); %Valeur du portefeuille de couverture actualisé
        V(i+1)=ValTheoriqueBS(t(i+1),S(i+1)); %Option calculée à l'aide de la formule de BS
    end
    
end


%On trace les graphs d'évolutions

function[]= Graph_Prtf_Couv()
    
    [Pactualise,V,A] = Simulation_Prtf_Couv();  
    
    %Graphe de l'évolution de l'option et du portefeuille de couverture
    %actualisé.
    figure;
    plot(t,Pactualise,t,V);
    title('Evolution de l option V et du portefeuille de couverture actualise Pactualise');
    xlabel('Temps t');
    ylabel('Prix');
    legend('Pactualise','V','Location','northeast');

    %Graphe de l'évolution des quantités d actions Ai et du cash Bi
    figure;
    plot(t,A,t,B);
    title('Courbes des quantités d actions Ai et du cash Bi');
    xlabel('Temps t');
    ylabel('Prix');
    legend('Ai : quantité d actions','Bi : Cash','Location','northwest');

    %Graphe d'erreur Pactualise - V
    figure;
    j=linspace(0,T,N+1);
    plot(j,Pactualise - V);
    title('Graphe d erreur de la couverture');
    xlabel('Temps t');
    ylabel('Erreur de couverture');
    
end

Graph_Prtf_Couv();


%% Calcul de la moyenne et de la variance de Profit&Loss final

%Definition de la valeur finale P&L

function[PandL]=ValeurFinale_PandL(Pactualise,V)
    PandL=Pactualise(N+1)-V(N+1);
end


%Pour chaque simulation du portefeuille de couverture, calcul de la valeur
%finale P&L

function[valeurFinale] = Calcul_PandL_valeurFinale(Nmc)
        for i = 1:Nmc
            [Pactualise,V] = Simulation_Prtf_Couv();
            valeurFinale(i)= ValeurFinale_PandL(Pactualise,V);
        end
end


%Espérance du P&L (moyenne arithmétique de P&L)

function[Esp_PandL]=Esperance_PandL()
    
    SommePandL=0;
    
    for k=1:Nmc
        [Pactualise,V,A]=Simulation_Prtf_Couv();
        SommePandL=SommePandL+ValeurFinale_PandL(Pactualise,V);     
    end
    
    Esp_PandL = SommePandL/Nmc;
    disp('Espérance du P&L = ');
    disp(Esp_PandL);
    
end

Esperance_PandL();


%Fonction de répartition du P&L

function[Fdr] = Fdr_PandL()
    
    X = Calcul_PandL_valeurFinale(Nmc);
    a=min(X);
    b=max(X);
    Nb = 100;
    x = linspace(a,b,Nb);
    
    for i=1:Nb
        compteur = 0;
        for j = 1:Nmc
            if (X(j)<=x(i))
                compteur = compteur + 1;
            end
            Fdr(i)=compteur/Nmc;
        end
    end
    
    figure;
    plot(x,Fdr,'b'); 
    title('Fonction de répartition de P&L');
    xlabel('Erreur de la couverture');
    ylabel('Fonction de répartition');
end

Fdr_PandL();


%Fonction de densité du P&L

function[densite] = Densite_PandL()
    
    X = Calcul_PandL_valeurFinale(Nmc);
    a=min(X);
    b=max(X);
    Nb = 100;
    deltax=(b-a)/Nb;
    x = linspace(a,b,Nb);
    
     for i = 1:Nb
         compteur = 0;
         for j=1:Nmc
            if (X(j)>x(i) && X(j)<=x(i+1))
                compteur = compteur + 1;
            end
            
         end
         proba(i)=compteur/Nmc;
         densite(i)=proba(i)/deltax/Nmc;
     end
     
     figure;
     plot(x,densite,'b');
     title('Densité de P&L');
     xlabel('Erreur de la couverture');
     ylabel('Fonction de densité');
end

Densite_PandL();



%% Rebalancement du portefeuille un nombre Ntrading de fois

%P&L final pour la quantité
%d'intervention du trading égale à Ntrading = 100,50,25,20,5,2,1
%ie. frequenceRebalancement = 1,2,4,5,20,50
%(pour Ntrading=1, l'option n'est pas couverte)


%------Simulation portefeuille et autres pour Ntrading différents-------%

%Rebalancement à chaque ti = Fonction Simulation_Prtf_Couv() Ligne 109


%Rebalancement tous les .... (frequence de rebalancement noté dans le code
%frequenceRebalancement)

function[Pactualise,V,A]=Simulation_Prtf_Couv_Ntrading(Ntrading)
    
    frequenceRebalancement = N/Ntrading;
    
    A(1)=Nx(d1(t(1),S(1)));
    P(1)=A(1)*S(1)+B(1);
    V(1)=ValTheoriqueBS(t(1),S(1));
    Pactualise(1)=V(1);
    
    modulo=frequenceRebalancement;
    
    for i=1:N
        
        if mod(i,modulo)==0
            S(i+1)=S(i)*exp((r-(Sigma^2)/2)*deltat+Sigma*sqrt(deltat)*randn); %Calcul du prix de l'action
            P(i+1)=A(i)*S(i+1)+(1+r*deltat)*B(i); %Valeur du portefeuille
            A(i+1)=Nx(d1(t(i+1),S(i+1))); %Quantité d'actions qu'il faut avoir pour couvrir une option
            B(i+1)=P(i+1)-A(i+1)*S(i+1); %Valeur du cash
            Pactualise(i+1)=P(i+1)+(V(1)-P(1))*exp(r*t(i+1)); %Valeur du portefeuille de couverture actualisé
            V(i+1)=ValTheoriqueBS(t(i+1),S(i+1)); %Option calculée à l'aide de la formule de BS
        else
            S(i+1)=S(i)*exp((r-(Sigma^2)/2)*deltat+Sigma*sqrt(deltat)*randn);
            P(i+1)=P(i);
            A(i+1)=A(i);
            B(i+1)=B(i);
            Pactualise(i+1)=Pactualise(i);
            V(i+1)=ValTheoriqueBS(t(i+1),S(i+1));
        end
        
    end
    
end


%-------Espérance du P&L (moyenne arithmétique de P&L)------------------%

%Rebalancement à chaque ti ie Ntrading=100, frequenceRebalancement = 1
function[Esp_PandL,Var_PandL]=Esp_Var_PandL()
    
    SommeEspPandL=0;
    SommeVarPandL=0;
    
    for k=1:Nmc
        [Pactualise,V,A]=Simulation_Prtf_Couv();
        SommeEspPandL=SommeEspPandL+ValeurFinale_PandL(Pactualise,V);
        SommeVarPandL=SommeVarPandL+(ValeurFinale_PandL(Pactualise,V))^2; 
    end
    
    Esp_PandL = SommeEspPandL/Nmc;
    Var_PandL = (SommeVarPandL/Nmc)-(Esp_PandL)^2;
    
    esp = ['Pour Ntrading = 100, Espérance du P&L = ',num2str(Esp_PandL)];
    var = ['Pour Ntrading = 100, Variance du P&L = ',num2str(Var_PandL)];
    disp('------------')
    disp(esp);
    disp(var);
    
end


%Rebalancement tous les.... (frequenceRebalancement)
function[Esp_PandL,Var_PandL]=Esp_Var_PandL_Ntrading(Ntrading)
    
    SommeEspPandL=0;
    SommeVarPandL=0;
    
    for k=1:Nmc
        [Pactualise,V,A]=Simulation_Prtf_Couv_Ntrading(Ntrading);
        SommeEspPandL=SommeEspPandL+ValeurFinale_PandL(Pactualise,V);
        SommeVarPandL=SommeVarPandL+(ValeurFinale_PandL(Pactualise,V))^2; 
    end
    
    Esp_PandL = SommeEspPandL/Nmc;
    Var_PandL = (SommeVarPandL/Nmc)-(Esp_PandL)^2;
    
    esp = ['Pour Ntrading = ',num2str(Ntrading), ', Espérance du P&L = ',num2str(Esp_PandL)];
    var = ['Pour Ntrading = ',num2str(Ntrading), ', Variance du P&L = ',num2str(Var_PandL)];
    disp('------------')
    disp(esp);
    disp(var);
    
end

Esp_Var_PandL(); %Ntrading=100
Esp_Var_PandL_Ntrading(50);
Esp_Var_PandL_Ntrading(25);
Esp_Var_PandL_Ntrading(20);
Esp_Var_PandL_Ntrading(5);
Esp_Var_PandL_Ntrading(2);


%----Tracer le graphe du nombre de parts (ratio de la couverture) Ai------%

function[]= Graph_Ratio_Ntrading(Ntrading)
    
    [Pactualise,V,A] = Simulation_Prtf_Couv_Ntrading(Ntrading);  
    
    figure;
    plot(t,A);
    title('Courbe du ratio de la couverture Ai');
    xlabel('Temps t');
    ylabel('Prix');
    
end

Graph_Ratio_Ntrading(10); %Ratio pour Ntrading=10 (Hedge 1 fois sur 10)


%----Tracer les graphes de fonctions de densité et de répartition---------%


%Pour chaque simulation du portefeuille de couverture, calcul de la valeur
%finale P&L selon Ntrading

function[valeurFinale] = Calcul_PandL_valeurFinale_Ntrading(Nmc,Ntrading)
        for i = 1:Nmc
            [Pactualise,V] = Simulation_Prtf_Couv_Ntrading(Ntrading);
            valeurFinale(i)= ValeurFinale_PandL(Pactualise,V);
        end
end


function[Fdr_PandL_1,Fdr_PandL_2,Fdr_PandL_4,Fdr_PandL_5,Fdr_PandL_20,Fdr_PandL_50] = Fdr_PandL_Ntrading()
    
    PandL1 = Calcul_PandL_valeurFinale(Nmc);
    PandL2 = Calcul_PandL_valeurFinale_Ntrading(Nmc,50);
    PandL4 = Calcul_PandL_valeurFinale_Ntrading(Nmc,25);
    PandL5 = Calcul_PandL_valeurFinale_Ntrading(Nmc,20);
    PandL20 = Calcul_PandL_valeurFinale_Ntrading(Nmc,5);
    PandL50 = Calcul_PandL_valeurFinale_Ntrading(Nmc,2);
    
    a1=min(PandL1);
    b1=max(PandL1);
    a2=min(PandL2);
    b2=max(PandL2);
    a4=min(PandL4);
    b4=max(PandL4);
    a5=min(PandL5);
    b5=max(PandL5);
    a20=min(PandL20);
    b20=max(PandL20);
    a50=min(PandL50);
    b50=max(PandL50);
    liste_a=[a1,a2,a4,a5,a20,a50];
    liste_b=[b1,b2,b4,b5,b20,b50];
    a=min(liste_a);
    b=max(liste_b);
    
    Nb = 100;
    x = linspace(a,b,Nb);
    
    for i=1:Nb
        
        compteur1 = 0;
        compteur2 = 0;
        compteur4 = 0;
        compteur5 = 0;
        compteur20 = 0;
        compteur50 = 0;
        
        for j = 1:Nmc
            
            if (PandL1(j)<=x(i))
                compteur1 = compteur1 + 1;
            end
            if (PandL2(j)<=x(i))
                compteur2 = compteur2 + 1;
            end
            if (PandL4(j)<=x(i))
                compteur4 = compteur4 + 1;
            end
            if (PandL5(j)<=x(i))
                compteur5 = compteur5 + 1;
            end
            if (PandL20(j)<=x(i))
                compteur20 = compteur20 + 1;
            end
            if (PandL50(j)<=x(i))
                compteur50 = compteur50 + 1;
            end
            
            Fdr_PandL_1(i)=compteur1/Nmc;
            Fdr_PandL_2(i)=compteur2/Nmc;
            Fdr_PandL_4(i)=compteur4/Nmc;
            Fdr_PandL_5(i)=compteur5/Nmc;
            Fdr_PandL_20(i)=compteur20/Nmc;
            Fdr_PandL_50(i)=compteur50/Nmc;
            
        end
    end
    
    figure;
    plot(x,Fdr_PandL_1,'b',x,Fdr_PandL_2,'y',x,Fdr_PandL_4,'r',x,Fdr_PandL_5,'g',x,Fdr_PandL_20,'k',x,Fdr_PandL_50,'m');
    title('Fonction de répartition P&L selon Ntrading');
    xlabel('Erreur de la couverture');
    ylabel('Fonction de répartition');
    legend('Rebalancement (Hedge) à chaque ti','Rebalancement 1 fois sur 2','Rebalancement 1 fois sur 4','Rebalancement 1 fois sur 5','Rebalancement 1 fois sur 20','Rebalancement 1 fois sur 50','Location','northwest');

end

Fdr_PandL_Ntrading();


function[Densite_PandL_1,Densite_PandL_2,Densite_PandL_4,Densite_PandL_5,Densite_PandL_20,Densite_PandL_50] = Densite_PandL_Ntrading()
    
    PandL1 = Calcul_PandL_valeurFinale(Nmc);
    PandL2 = Calcul_PandL_valeurFinale_Ntrading(Nmc,50);
    PandL4 = Calcul_PandL_valeurFinale_Ntrading(Nmc,25);
    PandL5 = Calcul_PandL_valeurFinale_Ntrading(Nmc,20);
    PandL20 = Calcul_PandL_valeurFinale_Ntrading(Nmc,5);
    PandL50 = Calcul_PandL_valeurFinale_Ntrading(Nmc,2);
    
    a1=min(PandL1);
    b1=max(PandL1);
    a2=min(PandL2);
    b2=max(PandL2);
    a4=min(PandL4);
    b4=max(PandL4);
    a5=min(PandL5);
    b5=max(PandL5);
    a20=min(PandL20);
    b20=max(PandL20);
    a50=min(PandL50);
    b50=max(PandL50);
    liste_a=[a1,a2,a4,a5,a20,a50];
    liste_b=[b1,b2,b4,b5,b20,b50];
    a=min(liste_a);
    b=max(liste_b);
    
    Nb = 100;
    deltax=(b-a)/Nb;
    x = linspace(a,b,Nb);
    
     for i = 1:Nb
         
        compteur1 = 0;
        compteur2 = 0;
        compteur4 = 0;
        compteur5 = 0;
        compteur20 = 0;
        compteur50 = 0;
         
         for j=1:Nmc
             
            if (PandL1(j)>x(i) && PandL1(j)<=x(i+1))
                compteur1 = compteur1 + 1;
            end

            if (PandL2(j)>x(i) && PandL2(j)<=x(i+1))
                compteur2 = compteur2 + 1;
            end
            if (PandL4(j)>x(i) && PandL4(j)<=x(i+1))
                compteur4 = compteur4 + 1;
            end

            if (PandL5(j)>x(i) && PandL5(j)<=x(i+1))
                compteur5 = compteur5 + 1; 
            end
            
            if (PandL20(j)>x(i) && PandL20(j)<=x(i+1))
                compteur20 = compteur20 + 1; 
            end
            
            if (PandL50(j)>x(i) && PandL50(j)<=x(i+1))
                compteur50 = compteur50 + 1; 
            end
            
         end
         
         proba1(i)=compteur1/Nmc;
         Densite_PandL_1(i)=proba1(i)/deltax/Nmc;
         proba2(i)=compteur2/Nmc;
         Densite_PandL_2(i)=proba2(i)/deltax/Nmc;
         proba4(i)=compteur4/Nmc;
         Densite_PandL_4(i)=proba4(i)/deltax/Nmc;
         proba5(i)=compteur5/Nmc;
         Densite_PandL_5(i)=proba5(i)/deltax/Nmc;
         proba20(i)=compteur20/Nmc;
         Densite_PandL_20(i)=proba20(i)/deltax/Nmc;
         proba50(i)=compteur50/Nmc;
         Densite_PandL_50(i)=proba50(i)/deltax/Nmc;
         
     end
     
     figure;
     plot(x,Densite_PandL_1,'b',x,Densite_PandL_2,'y',x,Densite_PandL_4,'r',x,Densite_PandL_5,'g',x,Densite_PandL_20,'k',x,Densite_PandL_50,'m');
     title('Fonction de densité de P&L selon Ntrading');
     xlabel('Erreur de la couverture');
     ylabel('Fonction de densité');
     legend('Rebalancement (Hedge) à chaque ti','Rebalancement 1 fois sur 2',' Rebalancement 1 fois sur 4','Rebalancement 1 fois sur 5','Rebalancement 1 fois sur 20','Rebalancement 1 fois sur 50','Location','northwest');

end

Densite_PandL_Ntrading();


%% Var : Value at Risk
%Voir la fonction de répartition, prendre y = 0.10 et déterminer x par
%report sur l'axe des abscisses


%VaR à 10% pour un hedging à chaque deltat (tracage de la fdr)
Fdr_PandL();


%VaR à 10% pour un hedging 1 fois sur 10 (tracage de la fdr)
%Je créer une fonction pour uniquement ce hedging en particulier

function[Fdr_PandL_10] = Fdr_PandL_Ntrading10()
    
    PandL10 = Calcul_PandL_valeurFinale_Ntrading(Nmc,10);

    a=min(PandL10);
    b=max(PandL10);
    
    Nb = 100;
    x = linspace(a,b,Nb);
    
    for i=1:Nb
        
        compteur10 = 0;
        
        for j = 1:Nmc
            
            if (PandL10(j)<=x(i))
                compteur10 = compteur10 + 1;
            end
            
            Fdr_PandL_10(i)=compteur10/Nmc;
            
        end
    end
    
    figure;
    plot(x,Fdr_PandL_10,'b');
    title('Fonction de répartition P&L - Hedging une fois sur 10');
    xlabel('Erreur de la couverture');
    ylabel('Fonction de répartition');

end

Fdr_PandL_Ntrading10();


end