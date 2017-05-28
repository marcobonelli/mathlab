
function [xBest, yBest, IGDbest, IGDmed, IGDworst] = alg_genetico(naval, optP, nobj, nexec)

% Aplicação aos Problemas DTLZ1 e DTLZ2 com 3 e 5 objetivos

% Entradas:
%  naval -> número de avaliações da função objetivo
%  optP  -> seleção do problema: 1 - DTLZ1 / Outros - DTLZ2
%  nobj  -> número de objetivos (3 ou 5 objetivos)
%  nexec -> número de execuções do algortimo para solução do problema

% Saídas:
%  xBest   -> matriz contendo as váriaveis dos individuos não dominados da execução com melhor IGD
%  yBest   -> matriz contendo a avaliação das funções objetivo para cada individuo de iBest
%  IGDbest -> melhor valor de IGD obtido (relativo a xBest)
%  IGDmed  -> média dos valores de IGD obtidos para as 'nexec' execuções
%  IGDwort -> pior valor de IGD obtido
    
%--------------------------------------------------------------------------

% calcula o número de individuos da população
if nobj == 3
    npop= 91;
else
    npop= 210;
end

% número de gerações
ngen= round(naval/npop);
%ngen

% número de variáveis
if optP == 1            % DTLZ1
    k = 5;    
else                    % DTLZ2
    k = 10;    
end

nvar = nobj + k - 1;

% restrições para as variáveis de busca 
% (delimita espaço de busca)
xmin= 0; xmax= 1;
ls = 0;

% Executa a solução do mesmo problema pelo número de vezes 'nexec' requerido: 
for Solucao=1:nexec
    %Solucao
    
    tic
    
      ge = 0;   % geração atual (inicial)

    % Passo 1: cria população inicial 
      P = rand(npop,nvar);
      
    %------------------------
    % Passo 2: avalia população inicial 
    if optP == 1   % DTLZ1
        
        P = dtlz1 (P,npop,nvar,nobj);
        
    else           % DTLZ2
        
        P = dtlz2 (P,npop,nvar,nobj);
    end
    
    % Ordena população por frentes de não-dominância 
    % & calcula distancia de multidão entre individuos da mesma frente
    P = FastNonDominatedSort(P,nobj,nvar,npop);
    P = CrowdingDistance(P,nobj,nvar,npop);
    

    % Passo 3: repete processo até atingir número máximo de gerações
     while ge <= ngen 

        ge = ge+1;          % atualiza geração atual    

        % Passo 3.1) Seleção
        Q = sgenetico(P, nvar, nobj, npop);

        % Passo 3.2) Cruzamento & Mutação
        Q = fgenetico(Q, xmin, xmax, npop, nvar, nobj);  

        % Passo 3.3) Avaliação 
        
        % Avalia população Qi
        if optP == 1   % DTLZ1

            Q = dtlz1 (Q,npop,nvar,nobj);

        else           % DTLZ2

            Q = dtlz2 (Q,npop,nvar,nobj);
        end

        % Passo 3.4) Junta populações de Pais e Filhos (elitismo)
            S = [P(:,1:(nvar+nobj));Q];     

        % Passo 3.5) Ordena individuos por frentes de não dominância
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
            
        % Passo 3.6) Calcula distância de multidão e seleciona 50% dos melhores indivíduos
        P = CrowdingDistance(S,nobj,nvar,npop);
                        
     end
    
    % verifica o número de soluções não dominadas finais obtidas
    if P(end,nvar+nobj+1) == 1
        nnd = npop;
    else
        nnd = find(P(:,nvar+nobj+1)>1,1)-1;   
    end
        
    sol_var = P(1:nnd,1:nvar);            % valores das variáveis de busca para a solução final
    sol_obj = P(1:nnd,nvar+1:nvar+nobj);  % variáveis da solução final
            
    % Calcula IGD das Soluções
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
    
    % a) melhor solução
    [IGDbest,id] = min(IGD);        
    xBest = sfinal(id).var;        % variáveis    
    yBest = sfinal(id).obj;        % objetivos 
    
    % IGD médio
    IGDmed = mean(IGD);
    
    % Pior IGD
    IGDworst = max(IGD);
end


%-----------------------------
% Funcao de Avaliação - DTLZ1
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
% Funcao de Avaliação - DTLZ2
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
% Função para avaliar a dominância entre duas soluções
function flag = domina(ObjSol1,ObjSol2)
      
% ObjSol1 -> vetor com valores de avaliação das funções objetos do indivíduo 1
% ObjSol2 -> vetor com valores de avaliação das funções objetos do indivíduo 2

    % Se as soluções são incomparáveis, retorna 0
    flag = 0;
    
    % Se houver relação de dominância, retorna o indice indicando o indivíduo dominante
    
    % a) solução 1 domina a solução 2
    if all(ObjSol2 >= ObjSol1) && any(ObjSol1 < ObjSol2)
       flag = 1;
        
    % b) solução 2 domina a solução 1
    elseif all(ObjSol1 >= ObjSol2) && any(ObjSol2 < ObjSol1)
       flag = 2;
    end
    
end

%-----------------------------
% Função para Ordenamento por Não Dominância
function xPoP = FastNonDominatedSort(xPoP, nobj, nvar, npop)

    % Indivíduos não dominados recebem rank = 1;
    % Indivíduos da segunda fronteira recebem rank = 2;
    % E assim sucessivamente...
    % Após ordenar fronteiras, calcula a distância de aglomeração para cada frente.
 
    % Non-Dominated Sort: Ordena a População Baseado em Não Dominância

    % inicializa primeira frente (não dominada)  
    frente = 1;
    F(frente).f = [];
    individuo = [];

    % posição para para guardar dados de rank (número da frente de dominância)
    pos = nobj + nvar + 1;      

    % 1) Compara a dominância entre todos os individuos da população, 
    %    dois a dois, e identifica a frente de não dominância.

    % para cada individuo "i" da população...
    for i = 1:npop

        individuo(i).n = 0;     % número de indivíduos que dominam "i" 
        individuo(i).p = [];    % conjunto que guardará todos indivíduos que "i" domina;

        % toma valores das funções objetivo para o indivíduo "i"
        Obj_Xi = xPoP(i,nvar+1:nvar+nobj);

        % para cada individuo "j" da população...
        for j = 1:npop

            % toma valores das funções objetivo para o indivíduo "j"
            Obj_Xj = xPoP(j,nvar+1:nvar+nobj);

            % verifica dominância
            flag = domina(Obj_Xi,Obj_Xj);

            % se "j" domina "i": incrementa o número de indivíduos que o dominam;
            if flag == 2
                individuo(i).n = individuo(i).n + 1;

            % então "i" domina "j": guarda índice do indivíduo "j" dominado por "i"
            elseif flag == 1
                individuo(i).p = [individuo(i).p j];
            end
        end   

        % se solução não for dominada por nenhuma outra... 
        % esta solução pertence a frente não dominada (rank=1)
        if individuo(i).n == 0
            xPoP(i, pos) = 1;                   % guarda rank na população
            F(frente).f = [F(frente).f i];      % salva individuo da frente não dominada
        end
    end

    % 2) Divide as demais soluções pelas frentes de dominância,
    %    conforme memória de dominação entre indivíduos da população

    % encontra as fronteiras seguintes:
    while ~isempty(F(frente).f)

       Qf = [];  % conjunto para guardar individuos da i-ésima fronteira;

       % para cada indivíduo "i" pertencente a última fronteira Fi de não dominância verificada...   
       for i = 1:length(F(frente).f)

           if ~isempty(individuo(F(frente).f(i)).p)

               % para cada um dos "j" indivíduos dominados pelos membros de Fi...
               for j = 1:length(individuo(F(frente).f(i)).p)

                   % decrementa o contador de dominação do indivíduo "j"
                   individuo(individuo(F(frente).f(i)).p(j)).n = individuo(individuo(F(frente).f(i)).p(j)).n - 1;

                   % verifica que nenhum dos individuos nas fronteiras subsequentes dominam "q"
                   if individuo(individuo(F(frente).f(i)).p(j)).n == 0

                        xPoP(individuo(F(frente).f(i)).p(j),pos) = frente + 1;   % guarda rank do indivíduo
                        Qf = [Qf individuo(F(frente).f(i)).p(j)];                % salva indivíduo da frente não dominada atual
                   end                
               end
           end
       end

       % atualiza posição para guardar indivíduos da próxima frente de dominância
       frente =  frente + 1;

       % salva individuos na frente de verificação corrente
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

% Função para Cálculo da Distância de Multidão
function xPoP = CrowdingDistance(xPoP, nobj, nvar, nsel)

    % Crowding Distance: Calcula distância de multidão 
    % (distância para as soluções vizinhas em cada dimensão do espaço de busca)

    % posição com informação do rank (número da frente de dominância)
    pos = nobj + nvar + 1;    

    % ordena individuos da população conforme nível de dominância
    [~,indice_fr] = sort(xPoP(:,pos));
    xPoP = xPoP(indice_fr,:);
       
    % verifica qual a ultima frente de dominância a entrar diretamente na população de pais
    lastF =  xPoP(nsel,pos);
    
    % verifica qual o pior rank de frente de não dominância na população
    worstF = max(xPoP(:,pos));
    
    % encontra a distância de multidão para cada indivíduo das frentes de dominância selecionadas
    for frente = 1:lastF

        % indice do primeiro individuo da fronteira
        if frente==1
            indice_ini = 1;    
        else
            indice_ini = indice_fim+1;
        end

        % indice do último indivíduo pertencente a fronteira
        if frente~=lastF || lastF < worstF
            indice_fim = find(xPoP(:,pos)>frente,1)-1;    
        else
            indice_fim = length(xPoP);
        end
        
        % número de soluções na frente
        nsolFi = indice_fim-indice_ini+1;

        % separa apenas as avaliações de função objetivo dos individuos na fronteira Fi
        xPop_Fi = xPoP(indice_ini:indice_fim, nvar+1:nvar+nobj);

        % inicializa vetor com valor nulo para as distâncias de multidão
        Di=zeros(1,nsolFi);

        % para cada objetivo "i"...
        for i = 1 : nobj

            % ordena indivíduos da fronteira baseado no valor do objetivo "i"  
            [~, indice_obj] = sort(xPop_Fi(:,i));
            
            % maior valor para o objetivo "i" - último indice
            f_max = xPop_Fi(indice_obj(end),i);

            % menor valor para o objetivo "i" - primeiro indice
            f_min = xPop_Fi(indice_obj(1),i);

            % atribui valor "infinito" para indivíduos na extremidade ótima da fronteira Fi
            Di(1,indice_obj(1)) = Di(1,indice_obj(1)) + Inf;

            % para individuos entre que não estão nas extremidades...
            for j = 2 : nsolFi

                % identifica valores da função objetivos das soluções vizinhas
                if j~=nsolFi
                    proximo   = xPop_Fi(indice_obj(j+1),i);
                else
                    % no extremo máximo, soma o semiperimetro apenas entre o único vizinho
                    proximo   = xPop_Fi(indice_obj(j),i);
                end

                anterior  = xPop_Fi(indice_obj(j-1),i);

                % calcula semi-perimetro normalizado
                Di(1,indice_obj(j)) = Di(1,indice_obj(j))+(proximo - anterior)/(f_max - f_min);
            end      
        end

        % guarda distâncias calculadas por individuo
        xPoP(indice_ini:indice_fim,pos+1)=Di';

        % ordena individuos da frente conforme distância de multidão
        [~,indice_fr] = sort(Di,'descend');
        xPoP(indice_ini:indice_fim,:) = xPoP(indice_fr+(indice_ini-1),:);
    end
    
    % seleciona apenas o valor 'nsel' de soluções
    xPoP = xPoP(1:nsel,:);
end

%--------------------

% Funcao para Cálculo da Distância Geracional Invertida (IGD)
function IGD = CalculaIGD(pareto, solucao)  
    % entradas:
    % pareto - Soluções da Fronteira de Pareto
    % solucao - Soluções não dominadas obtidas pelo algoritmo em teste

    % núm. de soluções da fronteira de Pareto
    [npareto,~] = size(pareto);
    
    % núm. de soluções obtidas pelo algoritmo desenvolvido
    [nsol,~] = size(solucao);
    
    dmin = zeros(1,npareto);    % guarda menores distâncias (di) entre a fronteira pareto e as soluções não dominadas
   
    d = zeros(npareto,nsol);    % dist. euclidiana entre cada ponto da fronteira pareto e cada solução não dominada
 
    % calcula distância euclidiana ponto a ponto
    for i=1:npareto
        for j=1:nsol            
            d(i,j) = norm(pareto(i,:)-solucao(j,:),2);
        end
        
        % guarda menor distância
        dmin(i) = min(d(i,:));
    end
    
    % realiza a média das menores distâncias
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
            % valor normalizado dos objetivos para uma solução "i" (aleatória)
            s = [xPop(i, 1) xPop(i, 2) xPop(i, 3)];

            % ponto de referência considerado
            w = [Zs(j) Zs(j) Zs(j)];

            % calculo da distância de projeção entre a solução e o ponto de referência
            dist(i, j) = norm(s-w*s'*w/norm(w)^2);
        end
    end
end
