%%TP3 - Vasicek Model
%%Part 2 - Simulation de différents cas de la courbe Yield

function[]=Partie_2()
    %% Fonctions définissant les valeurs de A, B, Yield et P
    function[B] = B_function(t,T,gamma)
        Tau = T-t;
        B = (1-exp(-gamma*Tau))/gamma;
    end


    function[A] = A_function(t,T,sigma,eta,gamma)
        Tau = T-t;
        A = (B_function(t,T,gamma)-Tau)*(eta*gamma-sigma^2/2)/gamma^2 - sigma^2*B_function(t,T,gamma)^2/4*gamma;
    end


    function[Y] = Yield(t, T, sigma, eta, gamma, r)
        Tau = T-t;
        Y = -(A_function(t,T,sigma,eta,gamma)-r*B_function(t,T,gamma))/Tau;
    end



    %% Fonctions pour la simulation du Modèle Vasicek selon différents cas et tracages des graphes


    %----------Cas 1 : Tracer Yield en fonction de r0 et T---------%

    %Valeurs
    t=0 ; 
    sigma = 0.02 ; 
    eta = 0.25*0.03 ;
    gamma = 0.25 ; 
    r0=linspace(0,1,100);
    T=linspace(1,30,30);

    function[]=SimulationCas1(t,T,sigma,eta,gamma,r0)

        %Initialisation
        Y=[];

        %Simulation du Yield selon T et r0
        for i = 1:30
            for j = 1:100
                Y(i,j) = Yield(t,T(i),sigma,eta,gamma,r0(j));
            end
        end

        %Tracé de la figure
        figure;
        surf(r0,T,Y)
        title("Graph de simulation Cas n°1 du modèle Vasicek selon différents paramètres")
        xlabel("r0")
        ylabel("T")
        zlabel("Yield")

    end

    %Appel de la fonction
    SimulationCas1(t,T,sigma,eta,gamma,r0)



    %----------Cas 2 : Tracer Yield en fonction de gamma et eta---------%

    %Valeurs
    t=0 ; 
    T = 10 ; 
    sigma = 0.02 ; 
    r0 = 0.1 ;
    eta = linspace(0,0.1,100) ;
    gamma = linspace(0.01,0.5,100) ; 

    function[]=SimulationCas2(t,T,sigma,eta,gamma,r0)

        %Initialisation
        Y=[];

        %Simulation du Yield selon eta et gamma
        for i = 1:100
            for j = 1:100
                Y(i,j) = Yield(t,T,sigma,eta(j),gamma(i),r0);
            end
        end

        %Tracé de la figure
        figure;
        surf(eta,gamma,Y)
        title("Graph de simulation Cas n°2 du modèle Vasicek selon différents paramètres")
        xlabel("eta")
        ylabel("gamma")
        zlabel("Yield")

    end

    %Appel de la fonction
    SimulationCas2(t,T,sigma,eta,gamma,r0)



    %----------Cas 3 : Tracer Yield en fonction de gamma et sigma---------%

    %Valeurs
    t=0 ; 
    T = 10 ; 
    r0 = 0.1 ;
    eta = 0.02 ;
    sigma = linspace(0,0.1,100) ; 
    gamma = linspace(0.01,0.5,100) ; 

    function[]=SimulationCas3(t,T,sigma,eta,gamma,r0)

        %Initialisation
        Y=[];

        %Simulation du Yield selon gamma et sigma
        for i = 1:100
            for j = 1:100
                Y(i,j) = Yield(t,T,sigma(j),eta,gamma(i),r0);
            end
        end

        %Tracé de la figure
        figure;
        surf(sigma,gamma,Y)
        title("Graph de simulation Cas n°3 du modèle Vasicek selon différents paramètres")
        xlabel("sigma")
        ylabel("gamma")
        zlabel("Yield")

    end

    %Appel de la fonction
    SimulationCas3(t,T,sigma,eta,gamma,r0)




    %---------Fonction simulant le Yield pour différentes valeurs de r-------%

    %Valeurs
    t=0 ;
    sigma = 0.02 ; 
    eta = 0.25*0.03 ;
    gamma = 0.25 ; 
    r1 = 0.01 ;
    r2 = 0.02 ;
    r3 = 0.035 ;
    r4 = 0.05 ;
    r0=[r1,r2,r3,r4] ;
    T = linspace(0,30,30) ;

    function[]=SimulationYield(t,T,sigma,eta,gamma,r0)

        %Initialisation
        Y=[];

        %Simulation du Yield selon 
        for j = 1:30
            Y1(j) = Yield(t,T(j),sigma,eta,gamma,r0(1));
            Y2(j) = Yield(t,T(j),sigma,eta,gamma,r0(2));
            Y3(j) = Yield(t,T(j),sigma,eta,gamma,r0(3));
            Y4(j) = Yield(t,T(j),sigma,eta,gamma,r0(4));
        end

        %Tracé de la figure
        figure;
        disp(Y1)
        plot(T,Y1,T,Y2,T,Y3,T,Y4)
        title("Simulation du Yield selon différentes valeurs de r0 et T")
        xlabel("T")
        ylabel("Yield")
        legend('r0 = 0.01','r0 = 0.02','r0 = 0.035', 'r0 = 0.05', 'Location','northeast');

    end

    %Appel de la fonction
    SimulationYield(t,T,sigma,eta,gamma,r0);


end
