function[]=Hedging_Transaction_Cost()
%%PARTIE PANDL ET PORTEFEUILLE AVEC COUT DE TRANSACTION

%% Valeurs des données fixes ou initiales

S(1)=100;
r=0.05;
Sigma=0.5; 
T=1;
N=260;
K=100;
t=linspace(0,T,N+1);
deltat=T/N; 
Nmc=1000;
k=0.001;
B(1)=100;

%%---------------------------------------------------------------------%%
%%-------------Hedging sans cout de transaction-------------------------%%
%%---------------------------------------------------------------------%%

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


%-------------------------Hedging a chaque deltat------------------------%


%           Simulation du portefeuille de couverture
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


%       Calcul de la moyenne et de la variance de Profit&Loss final

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




%%---------------------------------------------------------------------%%
%%-------------Hedging avec cout de transaction-------------------------%%
%%---------------------------------------------------------------------%%

function[f]=Nx_cout(x)
    f=0.5*(1+erf(x/sqrt(2)));
end

function[f]=d1_cout(t,S,Sigma,K)
    f=(log(S/K)+(r+(Sigma^2)/2)*(T-t))/(Sigma*sqrt(T-t));
end

function[f]=d2_cout(t,S,Sigma,K)
    f=(log(S/K)+(r-(Sigma^2)/2)*(T-t))/(Sigma*sqrt(T-t));
end

function[f]=ValTheoriqueBS_cout(t,S,Sigma,K)
    f=S*Nx_cout(d1_cout(t,S,Sigma,K))-K*exp(-r*(T-t))*Nx_cout(d2_cout(t,S,Sigma,K));
end


%           Simulation du portefeuille de couverture
function[P,V,A]=Simulation_Prtf_Couv_cout(deltat,K,k,Sigma)
    
    
    A(1)=Nx_cout(d1_cout(t(1),S(1),Sigma,K));
    V(1)=ValTheoriqueBS_cout(t(1),S(1),Sigma,K);
    P(1)=V(1) - k*abs(A(1))*S(1);
    B(1) = V(1) - A(1)*S(1) - k*abs(A(1))*S(1);

    for i=1:N
        
        S(i+1)=S(i)*exp((r-(Sigma^2)/2)*deltat+Sigma*sqrt(deltat)*randn); %Calcul du prix de l'action
        V(i+1)=ValTheoriqueBS_cout(t(i+1),S(i+1),Sigma,K); %Option calculée à l'aide de la formule de BS
        A(i+1)=Nx_cout(d1_cout(t(i+1),S(i+1),Sigma,K)); %Quantité d'actions qu'il faut avoir pour couvrir une option
        P(i+1)=A(i)*S(i+1)+B(i)*(1+r*deltat)-k*abs(A(i+1)-A(i))*S(i+1); %Valeur du portefeuille
        B(i+1)=B(i)*(1+r*deltat)+S(i+1)*(A(i)-A(i+1))-k*S(i+1)*abs(A(i)-A(i+1)); %Valeur du cash

    end
    
end


%On trace les graphs d'évolutions

function[]= Graph_Prtf_Couv_cout(deltat,K,k,Sigma)
    
    [P,V] = Simulation_Prtf_Couv_cout(deltat,K,k,Sigma);
    
    %Graphe de l'évolution de l'option et du portefeuille de couverture
    %actualisé.
    figure;
    plot(t,P,t,V);
    title('Evolution de l option V et du portefeuille de couverture P');
    xlabel('Temps t');
    ylabel('Prix');
    legend('P','V','Location','northeast');

end


%       Calcul de la moyenne et de la variance de Profit&Loss final


%Pour chaque simulation du portefeuille de couverture, calcul de la valeur
%finale P&L
function[valeurFinale] = Calcul_PandL_valeurFinale_cout(deltat,K,k,Sigma,Nmc)
        for i = 1:Nmc
            [P,V] = Simulation_Prtf_Couv_cout(deltat,K,k,Sigma);
            valeurFinale(i)= ValeurFinale_PandL(P,V);
        end
end


%Espérance du P&L (moyenne arithmétique de P&L)
function[Esp_PandL]=Esperance_PandL_cout(deltat,K,k,Sigma)
    
    SommePandL=0;
    
    for i=1:Nmc
        [P,V]=Simulation_Prtf_Couv_cout(deltat,K,k,Sigma);
        SommePandL=SommePandL+ValeurFinale_PandL(P,V);     
    end
    
    Esp_PandL = SommePandL/Nmc;
    disp('Espérance du P&L avec coût de transaction = ');
    disp(Esp_PandL);
    
end


%Fonction de répartition du P&L
function[Fdr] = Fdr_PandL_cout(deltat,K,k,Sigma,Nmc)
    
    X = Calcul_PandL_valeurFinale_cout(deltat,K,k,Sigma,Nmc);
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
    title('Fonction de répartition de P&L avec coût de transaction');
    xlabel('Erreur de la couverture');
    ylabel('Fonction de répartition');
end


%Fonction de densité du P&L
function[densite] = Densite_PandL_cout(deltat,K,k,Sigma,Nmc)
    
    X = Calcul_PandL_valeurFinale_cout(deltat,K,k,Sigma,Nmc);
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
     title('Densité de P&L avec coût de transaction');
     xlabel('Erreur de la couverture');
     ylabel('Fonction de densité');

end



%% TESTS : SIMULATIONS 

%1. delta t=1/260, rebalance every deltat, K=100, k=0.001
T=1;
N=260;
deltat=T/N;
t=linspace(0,T,N+1);
k=0.001;
K=100;
sigmaVol=Sigma*sqrt(1+(k/Sigma)*sqrt(2/(pi*deltat)));

Graph_Prtf_Couv_cout(deltat,K,k,sigmaVol);
Fdr_PandL_cout(deltat,K,k,sigmaVol,Nmc);
Densite_PandL_cout(deltat,K,k,sigmaVol,Nmc);
Fdr_PandL();
Densite_PandL();
Esperance_PandL_cout(deltat,K,k,sigmaVol)
Fdr_PandL_cout(deltat,K,k,sigmaVol,Nmc);


%2. delta t=1/1040, 1 fois sur 4, K=100, k=0.001
T=1;
N=1040;
deltat=T/N;
t=linspace(0,T,N+1);
k=0.001;
K=100;
sigmaVol=Sigma*sqrt(1+(k/Sigma)*sqrt(2/(pi*deltat)));

Graph_Prtf_Couv_cout(deltat,K,k,sigmaVol);
Fdr_PandL_cout(deltat,K,k,sigmaVol,Nmc);
Densite_PandL_cout(deltat,K,k,sigmaVol,Nmc);
Fdr_PandL();
Densite_PandL();
Esperance_PandL_cout(deltat,K,k,sigmaVol)
Fdr_PandL_cout(deltat,K,k,sigmaVol,Nmc);


%3. delta t=1/260, rebalance every deltat, K=100, k=0.01
T=1;
N=260;
deltat=T/N;
t=linspace(0,T,N+1);
k=0.01;
K=100;
sigmaVol=Sigma*sqrt(1+(k/Sigma)*sqrt(2/(pi*deltat)));

Graph_Prtf_Couv_cout(deltat,K,k,sigmaVol);
Fdr_PandL_cout(deltat,K,k,sigmaVol,Nmc);
Densite_PandL_cout(deltat,K,k,sigmaVol,Nmc);
Fdr_PandL();
Densite_PandL();
Esperance_PandL_cout(deltat,K,k,sigmaVol)
Fdr_PandL_cout(deltat,K,k,sigmaVol,Nmc);


%4. delta t=1/1040, rebalance 1 fois sur 4, K=120, k=0.01
T=1;
N=1040;
deltat=T/N;
t=linspace(0,T,N+1);
k=0.01;
K=120;
sigmaVol=Sigma*sqrt(1+(k/Sigma)*sqrt(2/(pi*deltat)));


Graph_Prtf_Couv_cout(deltat,K,k,sigmaVol);
Fdr_PandL_cout(deltat,K,k,sigmaVol,Nmc);
Densite_PandL_cout(deltat,K,k,sigmaVol,Nmc);
Fdr_PandL();
Densite_PandL();
Esperance_PandL_cout(deltat,K,k,sigmaVol)
Fdr_PandL_cout(deltat,K,k,sigmaVol,Nmc);



end