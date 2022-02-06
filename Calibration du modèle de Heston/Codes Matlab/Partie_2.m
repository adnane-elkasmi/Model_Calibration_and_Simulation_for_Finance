%TP5 - Calibration du modèle de Heston

%Partie 2 - Calibration du modèle de Heston

function[]= Partie_2()
    %% Valeurs initiales
    
    h1=0.1;
    h2=0.1;
    S0=10;
    eps=10^-4;
    lambda=0.01;
    r=0.01;
    k=3;
    rho=0.5;
    theta=0.2;
    eta=0.5;
    N=100;
    T=0.5;
    deltat=T/N;
    
    %Nombre de simulations Monte-Carlo
    Nmc=1000;
        
    %Simulations de lois Normales
    N1=randn(N,Nmc);
    N2=randn(N,Nmc);

    %Valeurs de marché initiales
    Kmarche=linspace(8,16,21);
    Vmarche=[2.0944,1.7488,1.4266,1.1456,0.8919,0.7068,0.5461,0.4187,0.3166,0.2425,0.1860,0.1370,0.0967,0.0715,0.0547,0.0381,0.0306,0.0239,0.0163,0.0139,0.0086];
    d=[1,1];

    
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


    
    function [Greek_Theta,Greek_Eta]=Calcul_Greek_Option(K,theta,eta,N1,N2)
        
        %Grecque "Theta"
        Prix_Option_Theta_Terme1=Prix_Call_Europeen_S0_fixe_Estimateur2(K,theta+h1,eta,N1,N2);
        Prix_Option_Theta_Terme2=Prix_Call_Europeen_S0_fixe_Estimateur2(K,theta-h1,eta,N1,N2);
        Greek_Theta=(Prix_Option_Theta_Terme1-Prix_Option_Theta_Terme2)/(2*h1);
        
        %Grecque "Eta"
        Prix_Option_Eta_Terme1=Prix_Call_Europeen_S0_fixe_Estimateur2(K,theta,eta+h2,N1,N2);
        Prix_Option_Eta_Terme2=Prix_Call_Europeen_S0_fixe_Estimateur2(K,theta,eta-h2,N1,N2);
        Greek_Eta=(Prix_Option_Eta_Terme1-Prix_Option_Eta_Terme2)/(2*h2);
        
    end
    
    
    
    %% Calibration du modèle Heston - Algorithme de Levenber-Marquardt
    
    function[]=Calibration_Heston()
        
        while norm(d)>eps 
            
            for p=1:21
                V(p)=Prix_Call_Europeen_S0_fixe_Estimateur2(Kmarche(p),theta,eta,N1,N2);

                Res(p)=Vmarche(p)-V(p);

                [GreekTheta,GreekEta]=Calcul_Greek_Option(Kmarche(p),theta,eta,N1,N2);

                J(p,1)=-GreekTheta;
                J(p,2)=-GreekEta;
            end
            
            Mat=J.'*J+lambda*eye(2);
            d=-inv(Mat)*J.'*Res.';
            theta = theta + d(1);
            eta = eta + d(2);

            if theta>1
                theta=0.2;
            end
            if theta<0
                theta=0.2;
            end

            if eta>1
                eta=0.5;
            end
            if eta<0
                eta=0.5;
            end

        end

        %Affichage des valeurs des paramètres estimés
        disp('Valeur du paramètre estimé theta = ');
        disp(theta);
        disp('Valeur du paramètre estimé eta = ');
        disp(eta)

        %Tracé de la Vmarche et VHeston en fonction de Kmarche
        figure;
        plot(Kmarche,Vmarche,'o',Kmarche,V);
        legend('V_{marche}','V_{heston}');
        title('Calibration du modèle Heston à T=T_{max}')
        xlabel('K_{marché}');
        ylabel('V');
        
    end

    %% Appel des fonctions / Tests
    
    Calibration_Heston();


end
