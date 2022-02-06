%TP5 - Calibration du modèle de Heston


%Partie 1 - Simulation du modèle de Heston
%Sous-Part 3 - Calcul des Grecques par Monte-Carlo


function[]=Partie_1_3()
    %% Valeurs initiales
    
    h1=0.1;
    h2=0.1;
    S0=10;
    T=0.5;
    N=100;
    r=0.1;
    k=3;
    rho=0.5;
    theta=0.2;
    eta=0.5;
    deltat=T/N;

    %Nombre de simulations Monte-Carlo
    Nmc=1000;
    
    %Simulations de lois Normales
    N1=randn(N,Nmc);
    N2=randn(N,Nmc);

    
    %% Implémentations des fonctions principales
    
    %Calcul de l'actif S et de la volatilité v
    function[S,v,S_sym,v_sym]=Calcul_ActifVola(theta,eta,N1,N2)

        S(1)=S0;
        v(1)=0.04;
        S_sym(1)=S0;
        v_sym(1)=0.04;

        for i=1:N
            v(i+1)=v(i)+k*(theta-v(i))*deltat+eta*sqrt(abs(v(i)))*sqrt(deltat)*N1(i)+(eta^2/4)*(deltat*N1(i)^2-deltat);
            S(i+1)=S(i)*exp((r-v(i)/2)*deltat+sqrt(abs(v(i)))*(rho*sqrt(deltat)*N1(i)+sqrt(1-rho^2)*sqrt(deltat)*N2(i)));
            v_sym(i+1)=v_sym(i)+k*(theta-v_sym(i))*deltat+eta*sqrt(abs(v_sym(i)))*sqrt(deltat)*(-N1(i))+(eta^2/4)*(deltat*(-N1(i))^2-deltat);
            S_sym(i+1)=S_sym(i)*exp((r-v_sym(i)/2)*deltat+sqrt(abs(v_sym(i)))*(rho*sqrt(deltat)*(-N1(i))+sqrt(1-rho^2)*sqrt(deltat)*(-N2(i))));
        end
        
    end


    %Calcul du pay off du call européen
    function[f] = Payoff_Call_Europeen(S,K)
        f = max(S-K, 0);
    end

    %Calcul du prix du call selon l'estimateur 2
    function[Prix] = Prix_Call_Europeen_S0_fixe_Estimateur2(K,theta,eta,N1,N2)
        
        Somme = 0;
        
        for n = 1 : Nmc
            N1_vec=N1(:,n);
            N2_vec=N2(:,n);
            [S,~,S_sym,~]=Calcul_ActifVola(theta,eta,N1_vec,N2_vec);
            Somme = Somme + Payoff_Call_Europeen(S(N+1),K) + Payoff_Call_Europeen(S_sym(N+1),K);
        end
        
        Prix = Somme/(2*Nmc);
        
    end


    %% Calcul des options Grecques Beta et Eta
    
    function []= Calcul_Greek_Option()
        
        %Calcul des prix des options Grecques
        for k = 1 : 41
            
            K(k) = (k-1) * 0.5;
            
            %Grecque "Theta"
            Prix_Option_Theta_Terme1(k)=Prix_Call_Europeen_S0_fixe_Estimateur2(K(k),theta+h1,eta,N1,N2);
            Prix_Option_Theta_Terme2(k)=Prix_Call_Europeen_S0_fixe_Estimateur2(K(k),theta-h1,eta,N1,N2);
            Greek_Theta(k)=(Prix_Option_Theta_Terme1(k)-Prix_Option_Theta_Terme2(k))/(2*h1);
            
            %Grecque "Eta"
            Prix_Option_Eta_Terme1(k)=Prix_Call_Europeen_S0_fixe_Estimateur2(K(k),theta,eta+h2,N1,N2);
            Prix_Option_Eta_Terme2(k)=Prix_Call_Europeen_S0_fixe_Estimateur2(K(k),theta,eta-h2,N1,N2);
            Greek_Eta(k)=(Prix_Option_Eta_Terme1(k)-Prix_Option_Eta_Terme2(k))/(2*h2);
   
        end

        %Affichage de la courbe du Grecque "Theta"
        figure;
        plot(K,Greek_Theta);
        xlabel('Valeur de K');
        ylabel('Valeur du grecque "Theta"');
        title('Courbe d evolution du grecque "Theta" en fonction de K');

        %Affichage de la courbe du Grecque "Theta"
        figure;
        plot(K,Greek_Eta);
        xlabel('Valeur de K');
        ylabel('Valeur du grecque "Eta"');
        title('Courbe d evolution du grecque "Eta" en fonction de K');

    end

    %% Appel des fonctions / Tests

    Calcul_Greek_Option();

end
