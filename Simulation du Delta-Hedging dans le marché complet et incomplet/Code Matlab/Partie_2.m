function[]=Partie_2()
%%PARTIE REPLIQUE OPTION

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


%% Simulation du portefeuille de couverture qui réplique l'option

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

%Calcul du portefeuille actualisé de départ :
V(1)=ValTheoriqueBS(t(1),S(1));
Pactualise(1)=V(1);
disp('Portefeuille actualisé de départ: P0actualise =  V0 =')
disp(Pactualise(1))

%Valeur du portefeuille de départ : 
P(1)=V(1);
disp('Portefeuille de départ: P0 = V0 =')
disp(P(1))

%Valeur du Cash
B(1)=V(1)-A(1)*S(1);
disp('Valeur du Cash du départ : B0 =')
disp(B(1))

%Une simulation du portefeuille de couverture et prix de l'option, ratio,
%etc...

function[Pactualise,V,A]=Simulation_Prtf_Couv_Option_Replique()
    
    A(1)=Nx(d1(t(1),S(1)));
    V(1)=ValTheoriqueBS(t(1),S(1));
    P(1)=V(1);
    B(1)=V(1)-A(1)*S(1);
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

function[]= Graph_Prtf_Couv_Option_Replique()
    
    [Pactualise,V,A] = Simulation_Prtf_Couv_Option_Replique();  
    
    %Graphe de l'évolution de l'option et du portefeuille de couverture
    %actualisé.
    figure;
    plot(t,Pactualise,t,V);
    title('Evolution de l option V et du portefeuille de couverture actualise Pactualise');
    xlabel('Temps t');
    ylabel('Prix');
    legend('Pactualise','V','Location','northeast');


    %Graphe d'erreur Pactualise - V
    figure;
    j=linspace(0,T,N+1);
    plot(j,Pactualise - V);
    title('Graphe d erreur de la couverture');
    xlabel('Temps t');
    ylabel('Erreur de couverture');
    
end

Graph_Prtf_Couv_Option_Replique();



end