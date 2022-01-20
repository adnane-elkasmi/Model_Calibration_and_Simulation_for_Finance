
%%Partie 2 - Calibration et reconstruction de la volatilite locale

function[] = Partie_2()
    %% Valeurs initiales définies

    %Valeurs utilisées
    Kmax = 20;
    S0 = 10;
    r = 0.1;
    Tmax = 0.5;
    N = 199;
    M = 49;
    epsilon = 10^-5;
    lambda = 0.001;
    d = [1,1];
    
    %Discrétisation de la valeur du strike K et de la maturité T
    K = linspace(0,Kmax,N+2);
    deltaK = Kmax/(N+1);
    t = linspace(0,Tmax,M+2);
    deltat = Tmax/(M+1);

    %Valeurs des paramètres des fonctions
    h = 0.01;
    beta1 = 1;
    beta2 = 1;
    beta = [beta1,beta2];
    
    Kmarche = linspace(7,14,15);
    Vmarche = [3.3634,2.9092,2.4703,2.0536,1.6666,1.3167,1.0100,0.7504,0.5389,0.3733,0.2491,0.1599,0.0986,0.0584,0.0332];

    %% Implémentation des fonctions 
    
    %-------------------------------------------------------------------%
    % Résolution de l'équation de Dupire et on détermine V
    % (fonction nécessaire pour la seconde fonction)
    %-------------------------------------------------------------------%
    
    function[V]=Prix_Dupire(h,beta)

        V=zeros(M+2,N+2);

        %Conditions initiales
        for i=1:N+2
            V(1,i)=max(S0-K(i),0);
            sigma(i)= h + beta(1)/(K(i)^beta(2));
        end

        %Conditions aux limites
        for n=1:M+2
            V(n,1)=S0;
            V(n,N+2)=0;
        end
        
        %Calcul de D, B et A
        for i=2:N+1
            D(i)=1 + (deltat/2)*((sigma(i)*K(i)/deltaK)^2) ;
            B(i)=-(deltat/4)*(K(i)*r/deltaK + ((sigma(i)*K(i))/deltaK)^2);
            A(i)=(deltat/4)*(K(i)*r/deltaK - ((sigma(i)*K(i))/deltaK)^2);
        end

        %Calcul de C
        for n=1:M+1
            for i=2:N+1
                if i==2
                    C(n,i)=-(deltat/4)*(K(i)*r/deltaK - (sigma(i)^2 * K(i)^2)/deltaK^2)*V(n,i+1) + (1-deltat/2*sigma(i)^2*K(i)^2/deltaK^2)*V(n,i) + deltat/4*(r*K(i)/deltaK + sigma(i)^2 * K(i)^2 /deltaK^2)*V(n,i-1) - S0*B(2);
                else
                    C(n,i)=-(deltat/4)*(K(i)*r/deltaK - (sigma(i)^2 * K(i)^2)/deltaK^2)*V(n,i+1) + (1-deltat/2*sigma(i)^2*K(i)^2/deltaK^2)*V(n,i) + deltat/4*(r*K(i)/deltaK + sigma(i)^2 * K(i)^2 /deltaK^2)*V(n,i-1);
                end
            end
            
            C_star(n,2)=C(n,2);
            D_star(2)=D(2);

            for i=3:N+1
                D_star(i)=D(i)-(B(i)*A(i-1))/D_star(i-1);
                C_star(n,i)=C(n,i)-(B(i)*C_star(n,i-1))/D_star(i-1);
            end

            %Calcul de V
            V(n+1,N+1)=C_star(n,N+1)/D_star(N+1);

            for i=N:-1:2
                V(n+1,i)=(C_star(n,i)-A(i)*V(n+1,i+1))/D_star(i);
            end      
            
        end
        
    end


    %-------------------------------------------------------------------%
    % Extraction des prix utiles parmi tout les prix calculés et stockés
    % dans V(n,i)
   %-------------------------------------------------------------------%
   
    function[V_DupireUtile,Vega_Utile,betaFinal,dFinal] = Prix_Dupire_Utiles(h,beta)
        
        while norm(d)>epsilon

        for p=1:15

            indice = Kmarche(p)/deltaK+1;

            V_Dupire = Prix_Dupire(h,beta);
            V_DupireUtile(p) = V_Dupire(51,indice);
            V2_Dupire = Prix_Dupire(0,beta);
            V2_DupireUtile(p) = V2_Dupire(51,indice);

            Vega_Utile(p) = (V_DupireUtile(p)-V2_DupireUtile(p))/h;

            Res(p) = Vmarche(p)-V_DupireUtile(p);
            J(p,1) = -Vega_Utile(p)/(Kmarche(p)^beta(2));
            J(p,2) = Vega_Utile(p)*log(Kmarche(p))*beta(1)/(Kmarche(p)^beta(2));
        end

        Mat = J.'*J+lambda*eye(2);
        d = -inv(Mat)*J.'*Res.';
        beta = beta + d.';

        end
        
        betaFinal = beta;
        dFinal = d;
        
    end



    %----------------------------------------------------------%
    %----Fonction permettant de calculer la Vega de l'option---%
    %----------------------------------------------------------%
    function[vega]=Vega_Dupire(h,beta)
        vega=(Prix_Dupire(h,beta)-Prix_Dupire(0,beta))/h;
    end




    %% Tests réalisés
    
    %Appel de la fonction
    [V_DupireUtile,Vega_Utile,betaFinal,dFinal] = Prix_Dupire_Utiles(h,beta);

    %Affichage de la valeur des paramètres Beta1 et Beta2 et de d après
    %calibration
    disp('La valeur de d est : ');
    disp(dFinal);
    disp(['La valeur du paramètre beta1 après calibration est : ', num2str(betaFinal(1))]);
    disp(['La valeur du paramètre beta2 après calibration est : ', num2str(betaFinal(2))]);
    
    
    
    
    %Visualisation de la surface de V avant calibration avec a=5, m=5
    V = Prix_Dupire(0,[1,1]);
    figure;
    mesh(K,t,V);
    title('Surface de V en fonction de T et K avant calibration');
    xlabel('K');
    ylabel('T');
    zlabel('V');

    %Visualisation de la surface de Vega avant calibration avec a=5, m=5
    Vega = Vega_Dupire(h,[1,1]);
    figure;
    mesh(K,t,Vega);
    title('Surface de Vega en fonction de T et K avant calibration');
    xlabel('K');
    ylabel('T');
    zlabel('Vega');
    
    %Visualisation de la surface de V après calibration
    V_calibrate = Prix_Dupire(0,betaFinal);
    figure;
    mesh(K,t,V_calibrate);
    title('Surface de V en fonction de T et K après calibration');
    xlabel('K');
    ylabel('T');
    zlabel('V calibré');

    %Visualisation de la surface de Vega après calibration
    Vega_calibrate = Vega_Dupire(h,betaFinal);
    figure;
    mesh(K,t,Vega_calibrate);
    title('Surface de Vega en fonction de T et K après calibration');
    xlabel('K');
    ylabel('T');
    zlabel('Vega calibré');
    

    
    
    %Visualisation de Vmarché et VDupire en 2D
    figure;
    plot(Kmarche,Vmarche,'b*',Kmarche,V_DupireUtile,'green');
    legend('V_{marche}','V_{Dupire}');
    title('Visualisation de Vmarché et VDupire en fonction de K')
    xlabel('K');
    ylabel('V');

    %Graphique de Vega en fonction de K 
    figure;
    plot(Kmarche,Vega_Utile,'b');
    title('Visualisation de Vega en fonction de K à maturité Tmax')
    xlabel('K');
    ylabel('Vega');

end