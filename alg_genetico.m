
function [xBest, yBest, IGDbest, IGDmed, IGDworst] = alg_genetico(naval, optP, nobj, nexec)

% Aplica��o aos Problemas DTLZ1 e DTLZ2 com 3 e 5 objetivos

% Entradas:
%  naval -> n�mero de avalia��es da fun��o objetivo
%  optP  -> sele��o do problema: 1 - DTLZ1 / Outros - DTLZ2
%  nobj  -> n�mero de objetivos (3 ou 5 objetivos)
%  nexec -> n�mero de execu��es do algortimo para solu��o do problema

% Sa�das:
%  xBest   -> matriz contendo as v�riaveis dos individuos n�o dominados da execu��o com melhor IGD
%  yBest   -> matriz contendo a avalia��o das fun��es objetivo para cada individuo de iBest
%  IGDbest -> melhor valor de IGD obtido (relativo a xBest)
%  IGDmed  -> m�dia dos valores de IGD obtidos para as 'nexec' execu��es
%  IGDwort -> pior valor de IGD obtido
    
%--------------------------------------------------------------------------

% calcula o n�mero de individuos da popula��o
if nobj == 3
    npop= 91;
else
    npop= 210;
end

% n�mero de gera��es
ngen= round(naval/npop);
%ngen

% n�mero de vari�veis
if optP == 1            % DTLZ1
    k = 5;    
else                    % DTLZ2
    k = 10;    
end

nvar = nobj + k - 1;

% restri��es para as vari�veis de busca 
% (delimita espa�o de busca)
xmin= 0; xmax= 1;
ls = 0;

% Executa a solu��o do mesmo problema pelo n�mero de vezes 'nexec' requerido: 
for Solucao=1:nexec
    %Solucao
    
    tic
    
      ge = 0;   % gera��o atual (inicial)

    % Passo 1: cria popula��o inicial 
      P = rand(npop,nvar);
      
    %------------------------
    % Passo 2: avalia popula��o inicial 
    if optP == 1   % DTLZ1
        
        P = dtlz1 (P,npop,nvar,nobj);
        
    else           % DTLZ2
        
        P = dtlz2 (P,npop,nvar,nobj);
    end
    
    % Ordena popula��o por frentes de n�o-domin�ncia 
    % & calcula distancia de multid�o entre individuos da mesma frente
    P = FastNonDominatedSort(P,nobj,nvar,npop);
    P = CrowdingDistance(P,nobj,nvar,npop);
    

    % Passo 3: repete processo at� atingir n�mero m�ximo de gera��es
     while ge <= ngen 

        ge = ge+1;          % atualiza gera��o atual    

        % Passo 3.1) Sele��o
        Q = sgenetico(P, nvar, nobj, npop);

        % Passo 3.2) Cruzamento & Muta��o
        Q = fgenetico(Q, xmin, xmax, npop, nvar, nobj);  

        % Passo 3.3) Avalia��o 
        
        % Avalia popula��o Qi
        if optP == 1   % DTLZ1

            Q = dtlz1 (Q,npop,nvar,nobj);

        else           % DTLZ2

            Q = dtlz2 (Q,npop,nvar,nobj);
        end

        % Passo 3.4) Junta popula��es de Pais e Filhos (elitismo)
            S = [P(:,1:(nvar+nobj));Q];     

        % Passo 3.5) Ordena individuos por frentes de n�o domin�ncia
            S = FastNonDominatedSort(S,nobj,nvar,npop*2);
        
         ls = ls + 1;
         if ls == 5
             sizeV = 0;
             for k = 1 : npop * 2
                 if S(k, nvar + nobj + 1) == 1
                     sizeV = sizeV + 1;
                     T(sizeV, :) = localsearch(S(k, 1:nvar), nvar, nobj, optP);
                 end
             end
             %T
             S = [S(:,1:(nvar+nobj));T]; 
             S = FastNonDominatedSort(S,nobj,nvar,npop*2+sizeV);
             %S
             ls = 0;
         end
            
        % Passo 3.6) Calcula dist�ncia de multid�o e seleciona 50% dos melhores indiv�duos
        P = CrowdingDistance(S,nobj,nvar,npop);
                        
     end
    
    % verifica o n�mero de solu��es n�o dominadas finais obtidas
    if P(end,nvar+nobj+1) == 1
        nnd = npop;
    else
        nnd = find(P(:,nvar+nobj+1)>1,1)-1;   
    end
        
    sol_var = P(1:nnd,1:nvar);            % valores das vari�veis de busca para a solu��o final
    sol_obj = P(1:nnd,nvar+1:nvar+nobj);  % vari�veis da solu��o final
            
    % Calcula IGD das Solu��es
    if optP == 1 && nobj==3        
       load('dtlz1_3d.mat');
       IGD(Solucao) = CalculaIGD(fronteiraReal, sol_obj);

    elseif optP == 1 && nobj==5        
       load('dtlz1_5d.mat');
       IGD(Solucao) = CalculaIGD(fronteiraReal, sol_obj);
       
    elseif optP ~= 1 && nobj==3          
       load('dtlz2_3d.mat');
       IGD(Solucao) = CalculaIGD(fronteiraReal, sol_obj);
       
    else
       load('dtlz2_5d.mat');
       IGD(Solucao) = CalculaIGD(fronteiraReal, sol_obj);       
    end   
        
    sfinal(Solucao).var = sol_var;
    sfinal(Solucao).obj = sol_obj;
    
end

    %  Retorna atributos requisitados:
    
    % a) melhor solu��o
    [IGDbest,id] = min(IGD);        
    xBest = sfinal(id).var;        % vari�veis    
    yBest = sfinal(id).obj;        % objetivos 
    
    % IGD m�dio
    IGDmed = mean(IGD);
    
    % Pior IGD
    IGDworst = max(IGD);
end


%-----------------------------
% Funcao de Avalia��o - DTLZ1
function xPop = dtlz1 (xPop,npop,nvar,nobj)
    
    k=5;
    f = zeros(1,nobj);
    
    for ind=1:npop

        x = xPop(ind,:);
    
        s = 0;

        for i = nobj:nvar
            s = s + (x(i)-0.5)^2 - cos(20*pi*(x(i)-0.5));
        end

        g = 100*(k+s);

        f(1) = 0.5 * prod(x(1:nobj-1)) * (1+g);
        
        for i = 2:nobj-1
            f(i) = 0.5 * prod(x(1:nobj-i)) * (1-x(nobj-i+1)) * (1+g);
        end

        f(nobj) = 0.5 * (1-x(1)) * (1+g);
        
        xPop(ind,nvar+1:nvar+nobj) = f;
        
    end
end

%-----------------------------
% Funcao de Avalia��o - DTLZ2
function xPop = dtlz2 (xPop,npop,nvar,nobj)
    
    k=10;
    f = zeros(1,nobj);

    for ind=1:npop
        
        x = xPop(ind,:);

        s = 0;
        
        for i = nobj:nvar
            s = s + (x(i)-0.5)^2;
        end
        
        g = s;

        cosx = cos(x*pi/2);
        sinx = sin(x*pi/2);

        f(1) =  (1+g) * prod(cosx(1:nobj-1));
        
        for i = 2:nobj-1
            f(i) = (1+g) * prod(cosx(1:nobj-i)) * sinx(nobj-i+1);
        end
        
        f(nobj) = (1+g) * sinx(1);
        xPop(ind,nvar+1:nvar+nobj) = f;
    end
end

%-----------------------------
% Fun��o para avaliar a domin�ncia entre duas solu��es
function flag = domina(ObjSol1,ObjSol2)
      
% ObjSol1 -> vetor com valores de avalia��o das fun��es objetos do indiv�duo 1
% ObjSol2 -> vetor com valores de avalia��o das fun��es objetos do indiv�duo 2

    % Se as solu��es s�o incompar�veis, retorna 0
    flag = 0;
    
    % Se houver rela��o de domin�ncia, retorna o indice indicando o indiv�duo dominante
    
    % a) solu��o 1 domina a solu��o 2
    if all(ObjSol2 >= ObjSol1) && any(ObjSol1 < ObjSol2)
       flag = 1;
        
    % b) solu��o 2 domina a solu��o 1
    elseif all(ObjSol1 >= ObjSol2) && any(ObjSol2 < ObjSol1)
       flag = 2;
    end
    
end

%-----------------------------
% Fun��o para Ordenamento por N�o Domin�ncia
function xPoP = FastNonDominatedSort(xPoP, nobj, nvar, npop)

    % Indiv�duos n�o dominados recebem rank = 1;
    % Indiv�duos da segunda fronteira recebem rank = 2;
    % E assim sucessivamente...
    % Ap�s ordenar fronteiras, calcula a dist�ncia de aglomera��o para cada frente.
 
    % Non-Dominated Sort: Ordena a Popula��o Baseado em N�o Domin�ncia

    % inicializa primeira frente (n�o dominada)  
    frente = 1;
    F(frente).f = [];
    individuo = [];

    % posi��o para para guardar dados de rank (n�mero da frente de domin�ncia)
    pos = nobj + nvar + 1;      

    % 1) Compara a domin�ncia entre todos os individuos da popula��o, 
    %    dois a dois, e identifica a frente de n�o domin�ncia.

    % para cada individuo "i" da popula��o...
    for i = 1:npop

        individuo(i).n = 0;     % n�mero de indiv�duos que dominam "i" 
        individuo(i).p = [];    % conjunto que guardar� todos indiv�duos que "i" domina;

        % toma valores das fun��es objetivo para o indiv�duo "i"
        Obj_Xi = xPoP(i,nvar+1:nvar+nobj);

        % para cada individuo "j" da popula��o...
        for j = 1:npop

            % toma valores das fun��es objetivo para o indiv�duo "j"
            Obj_Xj = xPoP(j,nvar+1:nvar+nobj);

            % verifica domin�ncia
            flag = domina(Obj_Xi,Obj_Xj);

            % se "j" domina "i": incrementa o n�mero de indiv�duos que o dominam;
            if flag == 2
                individuo(i).n = individuo(i).n + 1;

            % ent�o "i" domina "j": guarda �ndice do indiv�duo "j" dominado por "i"
            elseif flag == 1
                individuo(i).p = [individuo(i).p j];
            end
        end   

        % se solu��o n�o for dominada por nenhuma outra... 
        % esta solu��o pertence a frente n�o dominada (rank=1)
        if individuo(i).n == 0
            xPoP(i, pos) = 1;                   % guarda rank na popula��o
            F(frente).f = [F(frente).f i];      % salva individuo da frente n�o dominada
        end
    end

    % 2) Divide as demais solu��es pelas frentes de domin�ncia,
    %    conforme mem�ria de domina��o entre indiv�duos da popula��o

    % encontra as fronteiras seguintes:
    while ~isempty(F(frente).f)

       Qf = [];  % conjunto para guardar individuos da i-�sima fronteira;

       % para cada indiv�duo "i" pertencente a �ltima fronteira Fi de n�o domin�ncia verificada...   
       for i = 1:length(F(frente).f)

           if ~isempty(individuo(F(frente).f(i)).p)

               % para cada um dos "j" indiv�duos dominados pelos membros de Fi...
               for j = 1:length(individuo(F(frente).f(i)).p)

                   % decrementa o contador de domina��o do indiv�duo "j"
                   individuo(individuo(F(frente).f(i)).p(j)).n = individuo(individuo(F(frente).f(i)).p(j)).n - 1;

                   % verifica que nenhum dos individuos nas fronteiras subsequentes dominam "q"
                   if individuo(individuo(F(frente).f(i)).p(j)).n == 0

                        xPoP(individuo(F(frente).f(i)).p(j),pos) = frente + 1;   % guarda rank do indiv�duo
                        Qf = [Qf individuo(F(frente).f(i)).p(j)];                % salva indiv�duo da frente n�o dominada atual
                   end                
               end
           end
       end

       % atualiza posi��o para guardar indiv�duos da pr�xima frente de domin�ncia
       frente =  frente + 1;

       % salva individuos na frente de verifica��o corrente
       F(frente).f = Qf;

    end
end

function Q = sgenetico(xPop, nvar, nobj, npop)
    k = 0.75;

    Q = xPop;
    pnts1 = randperm(npop);
    pnts2 = randperm(npop);
    cnt = 1;
    pos = nobj + nvar + 1;
    
    for i = 1 : npop
        if xPop(pnts1(cnt), pos) < xPop(pnts2(cnt), pos)
            if rand() < k
                Q(cnt, :) = xPop(pnts1(cnt), :);
            else
                Q(cnt, :) = xPop(pnts2(cnt), :);
            end
        elseif xPop(pnts1(cnt), pos) > xPop(pnts2(cnt), pos)
            if rand() < k
                Q(cnt, :) = xPop(pnts2(cnt), :);
            else
                Q(cnt, :) = xPop(pnts1(cnt), :); 
            end
        else
            Q(cnt, :) = xPop(pnts1(cnt), :);
        end
        
        cnt = cnt + 1;
    end
end

%--------------------

% Fun��o para C�lculo da Dist�ncia de Multid�o
function xPoP = CrowdingDistance(xPoP, nobj, nvar, nsel)

    % Crowding Distance: Calcula dist�ncia de multid�o 
    % (dist�ncia para as solu��es vizinhas em cada dimens�o do espa�o de busca)

    % posi��o com informa��o do rank (n�mero da frente de domin�ncia)
    pos = nobj + nvar + 1;    

    % ordena individuos da popula��o conforme n�vel de domin�ncia
    [~,indice_fr] = sort(xPoP(:,pos));
    xPoP = xPoP(indice_fr,:);
       
    % verifica qual a ultima frente de domin�ncia a entrar diretamente na popula��o de pais
    lastF =  xPoP(nsel,pos);
    
    % verifica qual o pior rank de frente de n�o domin�ncia na popula��o
    worstF = max(xPoP(:,pos));
    
    % encontra a dist�ncia de multid�o para cada indiv�duo das frentes de domin�ncia selecionadas
    for frente = 1:lastF

        % indice do primeiro individuo da fronteira
        if frente==1
            indice_ini = 1;    
        else
            indice_ini = indice_fim+1;
        end

        % indice do �ltimo indiv�duo pertencente a fronteira
        if frente~=lastF || lastF < worstF
            indice_fim = find(xPoP(:,pos)>frente,1)-1;    
        else
            indice_fim = length(xPoP);
        end
        
        % n�mero de solu��es na frente
        nsolFi = indice_fim-indice_ini+1;

        % separa apenas as avalia��es de fun��o objetivo dos individuos na fronteira Fi
        xPop_Fi = xPoP(indice_ini:indice_fim, nvar+1:nvar+nobj);

        % inicializa vetor com valor nulo para as dist�ncias de multid�o
        Di=zeros(1,nsolFi);

        % para cada objetivo "i"...
        for i = 1 : nobj

            % ordena indiv�duos da fronteira baseado no valor do objetivo "i"  
            [~, indice_obj] = sort(xPop_Fi(:,i));
            
            % maior valor para o objetivo "i" - �ltimo indice
            f_max = xPop_Fi(indice_obj(end),i);

            % menor valor para o objetivo "i" - primeiro indice
            f_min = xPop_Fi(indice_obj(1),i);

            % atribui valor "infinito" para indiv�duos na extremidade �tima da fronteira Fi
            Di(1,indice_obj(1)) = Di(1,indice_obj(1)) + Inf;

            % para individuos entre que n�o est�o nas extremidades...
            for j = 2 : nsolFi

                % identifica valores da fun��o objetivos das solu��es vizinhas
                if j~=nsolFi
                    proximo   = xPop_Fi(indice_obj(j+1),i);
                else
                    % no extremo m�ximo, soma o semiperimetro apenas entre o �nico vizinho
                    proximo   = xPop_Fi(indice_obj(j),i);
                end

                anterior  = xPop_Fi(indice_obj(j-1),i);

                % calcula semi-perimetro normalizado
                Di(1,indice_obj(j)) = Di(1,indice_obj(j))+(proximo - anterior)/(f_max - f_min);
            end      
        end

        % guarda dist�ncias calculadas por individuo
        xPoP(indice_ini:indice_fim,pos+1)=Di';

        % ordena individuos da frente conforme dist�ncia de multid�o
        [~,indice_fr] = sort(Di,'descend');
        xPoP(indice_ini:indice_fim,:) = xPoP(indice_fr+(indice_ini-1),:);
    end
    
    % seleciona apenas o valor 'nsel' de solu��es
    xPoP = xPoP(1:nsel,:);
end

%--------------------

% Funcao para C�lculo da Dist�ncia Geracional Invertida (IGD)
function IGD = CalculaIGD(pareto, solucao)  
    % entradas:
    % pareto - Solu��es da Fronteira de Pareto
    % solucao - Solu��es n�o dominadas obtidas pelo algoritmo em teste

    % n�m. de solu��es da fronteira de Pareto
    [npareto,~] = size(pareto);
    
    % n�m. de solu��es obtidas pelo algoritmo desenvolvido
    [nsol,~] = size(solucao);
    
    dmin = zeros(1,npareto);    % guarda menores dist�ncias (di) entre a fronteira pareto e as solu��es n�o dominadas
   
    d = zeros(npareto,nsol);    % dist. euclidiana entre cada ponto da fronteira pareto e cada solu��o n�o dominada
 
    % calcula dist�ncia euclidiana ponto a ponto
    for i=1:npareto
        for j=1:nsol            
            d(i,j) = norm(pareto(i,:)-solucao(j,:),2);
        end
        
        % guarda menor dist�ncia
        dmin(i) = min(d(i,:));
    end
    
    % realiza a m�dia das menores dist�ncias
    IGD = mean(dmin);
end

function Q = fgenetico(xPop, xmin, xmax, npop, nvar, nobj)
    Q = zeros(npop,nvar);
    [Zs,nZs] = ptos_referencia(nobj);
    dist = p_ortogonal(xPop, Zs, nZs, npop);
    dref = zeros(npop, 2);
    
    for k = 1 : nZs
        dref(:, 1) = 1 : npop;
        dref(:, 2) = dist(:, k);
        [~, indice_fr] = sort(dist(:, 2));
        x2Pop = xPop(indice_fr, :);
        nref = ceil(npop * 0.05);
        for i = 1 : nvar
            xm = mean(x2Pop(1:nref, i));
            xd = std(x2Pop(1:nref, i));
            var = normrnd(xm, xd);
            while var < xmin || var > xmax
                var = normrnd(xm, xd);
            end
            Q(k, i) = var;
        end
    end
end

function pPop = localsearch(p, nvar, nobj, optP)
    xmax = 1;
    xmin = 0;

    if optP == 1   % DTLZ1
        p = dtlz1 (p, 1, nvar, nobj);
    else           % DTLZ2
        p = dtlz2 (p, 1, nvar, nobj);
    end
    
    pfo = 0;
    for i = 1 : nobj
        pfo = pfo + (1/nobj)*p(nvar + i);
    end
    
    qPop = zeros(9, nvar);
    i = 1;
    while i < nvar + 1
        dxmax = xmax - p(i);
        dxmin = p(i) - xmin;
        
        for k = 1 : 4
           qPop(k, i) = p(i) - dxmin * (k / 4);
        end
        qPop(5, i) = p(i);
        for k = 1 : 4
           qPop(5 + k, i) = p(i) + dxmax * (k / 4);
        end
        
        for k = 1 : nvar
            if k ~= i
                qPop(:, k) = p(k);
            end
        end
        
        if optP == 1   % DTLZ1
            qPop = dtlz1 (qPop, 9, nvar, nobj);
        else           % DTLZ2
            qPop = dtlz2 (qPop, 9, nvar, nobj);
        end
        
        qfo = zeros(9, 1);
        for j = 1 : nobj
            for z = 1 : 9
                qfo(z, 1) = qfo(z, 1) + (1/nobj)*qPop(z, nvar + j);
            end
        end
        
        if min(qfo(:, 1)) >= pfo
            i = i + 1;
        else
            [~, indice_fr] = sort(qfo(:, 1));
            p = qPop(indice_fr(1), :);
            i = 1;
            pfo = min(qfo(:, 1));
        end
    end
    
    pPop = p;
end

function dist = p_ortogonal(xPop, Zs, nZs, npop)
    dist = zeros(npop, nZs);

    for i = 1 : npop
        for j = 1 : nZs
            % valor normalizado dos objetivos para uma solu��o "i" (aleat�ria)
            s = [xPop(i, 1) xPop(i, 2) xPop(i, 3)];

            % ponto de refer�ncia considerado
            w = [Zs(j) Zs(j) Zs(j)];

            % calculo da dist�ncia de proje��o entre a solu��o e o ponto de refer�ncia
            dist(i, j) = norm(s-w*s'*w/norm(w)^2);
        end
    end
end
