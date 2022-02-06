%TP5 - Calibration du modèle de Heston


%Partie 1 - Simulation du modèle de Heston
%Sous-Part 2 - Calcul du prix du Call par Monte Carlo


function[]=Partie_1_2()
    %% Valeurs initiale
    
    K=10; %Pour que la courbe soit centrée en K = 10
    r=0.01;
    k=2;
    rho=0;
    theta=0.04;
    eta=0.3;
    T=0.5;
    N=100;
    deltat=T/N;
    
    %Nombre de simulations Monte-Carlo
    Nmc=1000;
    

    %% Fonctions initiales à implémenter
    
    %Calcul de l'actif S et de la volatilité v
    function[S,v,S_sym,v_sym]=Calcul_ActifVola(S0)
        
        S(1)=S0;
        v(1)=0.04;
        S_sym(1)=S0;
        v_sym(1)=0.04;
        
        for i=1:N
            N1=randn;
            N2=randn;
            v(i+1)=v(i)+k*(theta-v(i))*deltat+eta*sqrt(abs(v(i)))*sqrt(deltat)*N1+(eta^2/4)*(deltat*N1^2-deltat);
            S(i+1)=S(i)*exp((r-v(i)/2)*deltat+sqrt(abs(v(i)))*(rho*sqrt(deltat)*N1+sqrt(1-rho^2)*sqrt(deltat)*N2));
            v_sym(i+1)=v_sym(i)+k*(theta-v_sym(i))*deltat+eta*sqrt(abs(v_sym(i)))*sqrt(deltat)*(-N1)+(eta^2/4)*(deltat*(-N1)^2-deltat);
            S_sym(i+1)=S_sym(i)*exp((r-v_sym(i)/2)*deltat+sqrt(abs(v_sym(i)))*(rho*sqrt(deltat)*(-N1)+sqrt(1-rho^2)*sqrt(deltat)*(-N2)));
        end
        
    end

    %Calcul du pay off du call européen
    function[f] = Payoff_Call_Europeen(S)
        f = max(S-K, 0);
    end


    %Calcul du prix du call selon l'estimateur utilisé
    function[Prix1, Prix2] = Prix_Call_Europeen_S0_fixe(S0)
        
        Somme1 = 0;
        Somme2 = 0;
        
        for n = 1 : Nmc
            [S,~,S_sym,~]=Calcul_ActifVola(S0);
            Somme1 = Somme1 + Payoff_Call_Europeen(S(N+1));
            Somme2 = Somme2 + Payoff_Call_Europeen(S(N+1)) + Payoff_Call_Europeen(S_sym(N+1));
        end
        
        Prix1 = Somme1/Nmc; %Estimateur 1
        Prix2 = Somme2/(2*Nmc); %Estimateur 2
        
    end


    %% Fonction pour affichage de la courbe du call selon les deux estimateurs
    
    function [Prix_estim_1, Prix_estim_2]= Courbe_Call_Europeen()
        
        for k = 1 : 200
            S0(k) = (k-1) * 0.1;
            [Prix_estim_1(k), Prix_estim_2(k)] = Prix_Call_Europeen_S0_fixe(S0(k));
        end

        figure;
        plot(S0,Prix_estim_1,S0,Prix_estim_2,S0,Payoff_Call_Europeen(S0));
        legend('Prix du call avec estimateur 1','Prix du call avec estimateur 2','Pay-off du call');
        xlabel('Valeur de S0');
        ylabel('Prix de l option');
        title('Prix du Call Europeen en fonction de S0');

    end


    %% Appel des fonctions / Tests
    
    %Prix de du call pour K = 10
    disp('Prix du Call Européen avec Monte-Carlo pour K = 10 : ')
    disp(Prix_Call_Europeen_S0_fixe(10));

    %Tracé des courbes du call selon l'estimateur utilisé
    Courbe_Call_Europeen();
    
end
