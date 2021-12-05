function[]=Partie_3()
%Partie GAMMA HEDGING

%% Valeurs des données fixes ou initiales

S(1)=1;
r=0.05;
Sigma=0.5;
T1=5; %Maturité du call à couvrir
T2=10; %Maturité du call de l'instrument de couverture
N=100;
t=linspace(0,T1,N+1);
deltat=T1/N; 
Nmc=1000;
K1=1.5;
K2=1.5;


%% Simulation de la couverture d'une option V(S,t) à l'aide des options du type C(S,t), des actions Ai et du cash Bi

%Option V (call a couvrir) --> indice 1
%Option C (call de l'instrument de couverture) --> indice 2

%Implémentations de toutes les fonctions nécessaires

function[f]=GammaV(t,S)
    f=exp(-0.5*(d11(t,S))^2)/(S*Sigma*sqrt(2*pi*(T1-t)));  
end

function[f]=GammaC(t,S)
    f=exp(-0.5*(d12(t,S))^2)/(S*Sigma*sqrt(2*pi*(T2-t)));  
end

function[f]=Nx(x)
    f=0.5*(1+erf(x/sqrt(2)));
end

function[f]=d11(t,S)
    f=(log(S/K1)+(r+(Sigma^2)/2)*(T1-t))/(Sigma*sqrt(T1-t));
end

function[f]=d12(t,S)
    f=(log(S/K2)+(r+(Sigma^2)/2)*(T2-t))/(Sigma*sqrt(T2-t));
end

function[f]=d21(t,S)
    f=(log(S/K1)+(r-(Sigma^2)/2)*(T1-t))/(Sigma*sqrt(T1-t));
end

function[f]=d22(t,S)
    f=(log(S/K2)+(r-(Sigma^2)/2)*(T2-t))/(Sigma*sqrt(T2-t));
end

function[f]=ValTheoriqueBS(t,S)
    f=S*Nx(d11(t,S))-K1*exp(-r*(T1-t))*Nx(d21(t,S));
end

function[f]=C_BS(t,S)
    f=S*Nx(d12(t,S))-K2*exp(-r*(T2-t))*Nx(d22(t,S));
end


%----------------------------------------------------------------%
%Calcul de la quantité initiale des options du type C
G(1)=GammaV(t(1),S(1))/GammaC(t(1),S(1));
disp('Quantité initiale des options du type C : G0 = ') 
disp(G(1))

%Calcul de la quantité initiale des actions
disp('Quantité initiale des actions : A0 = ') 
A(1)=Nx(d11(t(1),S(1)))-G(1)*Nx(d12(t(1),S(1)));
disp(A(1))

%Valeur du Cash
B(1)=1;

%Valeur du portefeuille de départ : 
P(1)=A(1)*S(1)+G(1)*C_BS(t(1),S(1))+B(1);
disp('Portefeuille de départ: P0 = ')
disp(P(1))

%Calcul du portefeuille actualisé de départ :
V(1)=ValTheoriqueBS(t(1),S(1));
Pactualise(1)=V(1); 
disp('Portefeuille actualisé de départ: P0actualise =  V0 =')
disp(Pactualise(1))


%Simulation du portefeuille de couverture - Gamma Hedging

function[Pactualise,V,A,C]=Simulation_Prtf_Couv_Gamma_Hedging()
    
    %Valeurs initiales
    
    G(1)=GammaV(t(1),S(1))/GammaC(t(1),S(1));
    A(1)=Nx(d11(t(1),S(1)))-G(1)*Nx(d12(t(1),S(1)));
    B(1)=1;
    V(1)=ValTheoriqueBS(t(1),S(1));
    C(1)=C_BS(t(1),S(1));
    P(1)=A(1)*S(1)+G(1)*C(1)+B(1);
    Pactualise(1)=V(1);
    
    %Recalculs et simulations
    
    for i=1:N
        
        S(i+1)=S(i)*exp((r-(Sigma^2)/2)*deltat+Sigma*sqrt(deltat)*randn); %Calcul du prix de l'action
        P(i+1)=A(i)*S(i+1)+G(i)*C_BS(t(i+1),S(i+1))+(1+r*deltat)*B(i); %Valeur du portefeuille
        G(i+1)=GammaV(t(i+1),S(i+1))/GammaC(t(i+1),S(i+1)); %Quantité d'options qu'il faut pour couvrir l'option V
        A(i+1)=Nx(d11(t(i+1),S(i+1)))-G(i+1)*Nx(d12(t(i+1),S(i+1))); %Quantité d'actions qu'il faut avoir pour couvrir une option
        B(i+1)=P(i+1)-A(i+1)*S(i+1)-G(i+1)*C_BS(t(i+1),S(i+1)); %Valeur du cash
        Pactualise(i+1)=P(i+1)+(V(1)-P(1))*exp(r*t(i+1)); %Valeur du portefeuille de couverture actualisé
        V(i+1)=ValTheoriqueBS(t(i+1),S(i+1)); %Option V calculée à l'aide de la formule de BS
        C(i+1)=C_BS(t(i+1),S(i+1)); %Option C calculée à l'aide de la formule de BS
        
    end
    
end


%On trace les graphs d'évolutions

function[]= Graph_Prtf_Couv_Gamma_Hedging()
    
    [Pactualise,V,A,C] = Simulation_Prtf_Couv_Gamma_Hedging();  
    
    %Graphe de l'évolution de l'option c et du portefeuille de couverture
    %actualisé.
    figure;
    plot(t,Pactualise,t,C);
    title('Evolution de l option C et du portefeuille de couverture actualise Pactualise');
    xlabel('Temps t');
    ylabel('Prix');
    legend('Pactualise','C','Location','northeast');
    
    %Graphe de l'évolution de l'option V et du portefeuille de couverture
    %actualisé.
    figure;
    plot(t,Pactualise,t,V);
    title('Evolution de l option V et du portefeuille de couverture actualise Pactualise');
    xlabel('Temps t');
    ylabel('Prix');
    legend('Pactualise','V','Location','northeast');

    %Graphe d'erreur Pactualise - V
    figure;
    j=linspace(0,T1,N+1);
    plot(j,Pactualise - V);
    title('Graphe d erreur de la couverture');
    xlabel('Temps t');
    ylabel('Erreur de couverture');
    
end

Graph_Prtf_Couv_Gamma_Hedging();



%% Calcul de la moyenne et de la variance du P&L final


%Definition de la valeur finale P&L

function[PandL]=ValeurFinale_PandL_Gamma_Hedging(Pactualise,V)
    PandL=Pactualise(N+1)-V(N+1);
end


%Pour chaque simulation du portefeuille de couverture, calcul de la valeur
%finale P&L

function[valeurFinale] = Calcul_PandL_valeurFinale_Gamma_Hedging(Nmc)
        for i = 1:Nmc
            [Pactualise,V] = Simulation_Prtf_Couv_Gamma_Hedging();
            valeurFinale(i)= ValeurFinale_PandL_Gamma_Hedging(Pactualise,V);
        end
end


%Espérance du P&L (moyenne arithmétique de P&L)

function[Esp_PandL]=Esperance_PandL_Gamma_Hedging()
    
    SommePandL=0;
    
    for k=1:Nmc
        [Pactualise,V,A]=Simulation_Prtf_Couv_Gamma_Hedging();
        SommePandL=SommePandL+ValeurFinale_PandL_Gamma_Hedging(Pactualise,V);     
    end
    
    Esp_PandL = SommePandL/Nmc;
    disp('Espérance du P&L = ');
    disp(Esp_PandL);
    
end

Esperance_PandL_Gamma_Hedging();


%Fonction de répartition du P&L

function[Fdr] = Fdr_PandL_Gamma_Hedging()
    
    X = Calcul_PandL_valeurFinale_Gamma_Hedging(Nmc);
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
    title('Fonction de répartition de P&L - Gamma Hedging');
    xlabel('Erreur de la couverture');
    ylabel('Fonction de répartition');
end

Fdr_PandL_Gamma_Hedging();


%Fonction de densité du P&L

function[densite] = Densite_PandL_Gamma_Hedging()
    
    X = Calcul_PandL_valeurFinale_Gamma_Hedging(Nmc);
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
     title('Densité de P&L - Gamma Hedging');
     xlabel('Erreur de la couverture');
     ylabel('Fonction de densité');
end

Densite_PandL_Gamma_Hedging();



end
