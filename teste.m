    %   Script de projeto final da disciplina de Projeto de Filtros Digitais, o
%   qual é responsável por testar cantos de pássaros, e assim compará-lo com
%   o banco de dados. 
%
%   sample_teste: áudios a serem testados com cantos de pássaros
%   sample_treino: áudios da base de dados de cada canto (armazenado como
%   DBx.wav)
%   sample_imagens: imagens das espécies conforme sample_treino
%   savesteste: salva os arquivos desse script para comparação com database
%
%   Data: 12/10/19
%   Nota: Manter o database atualizado.
%
%   Autor: João Brum
%   

clc
clear all
close all

warning ('off');

% Seleciona o arquivo na pasta para testar o canto
[file,path] = uigetfile('*.wav');
 
% Armazena o áudio para treinamento (5 segundos equivalente em amostras)
[y,Fs] = audioread(fullfile(path,file));

%Reproduz 5 segundos do audio a ser testado
sound(y(1:220500,1),Fs);

pause(4);

% Retira o nível DC
y1 = y(:,1)-mean(y(:,1));

%% Decimador
% Fator de decimação
DecFactor = 3;
h = fir1(30,1/DecFactor,'low');
Hm = mfilt.firdecim(DecFactor,h);
Y1 = filter(Hm,y1);

% Frequência de amostragem obtida pela decimação
Fsnew = Fs/DecFactor;

% Vetor de tempo (amostras->segundos)
t = 0:1/Fsnew:(length(Y1)-1)*(1/Fsnew);
% Plot do dominio tempo depois da decimação
% figure(1);
% plot(t,Y1);
% title 'Domínio tempo pós decimação'
% xlabel 'Tempo(s)'
% ylabel 'Amplitude (V)'

% Plot do dominio frequência 
% figure(2);
% plotfft(Y1,Fsnew);
% title 'Domínio frequência pós decimação'
% Divide o sinal no domínio tempo em 40ms
frame = round(0.04*Fsnew);

% Número de segmentos de 40ms
segmentos = floor(length(Y1)/frame);

% Número de bandas
n_filtros=30;
% Filtra o sinal conforme os 30 filtros.
[freqSignals,n] = filterbank(Y1,Fsnew,t,n_filtros);

% Inicializa os parâmetros
yy1 = zeros(segmentos,n);
yy2 = zeros(segmentos,n);

  for j = 1:n
%              figure(2)
%              title(sprintf('Frequência Dominante em cada segmento - Banda %d',j));
%              xlabel('Segmento nº');
%              ylabel('Frequência');
%              hold on;
       for i = 1:segmentos
            % Segmento em questão
            A2 = freqSignals(((i-1)*frame+1:i*frame),j);
            % FFT do segmento em questão
            F2 = abs(fft(A2));
            F2 = F2(1:length(F2)/2);
            F2 = F2/length(F2);
            % Vetor de frequências correspondente
            f2 = linspace(0,pi*(Fsnew/(2*pi)),length(F2));
            % Componente e respectiva frequência dominante
            [val2,ind2] = max(F2);
            Freqs2 = f2(ind2);
            % Componente da frequência dominante
            yy1(i,j) =  val2;
            % energia média
            M1 = mean(abs(A2.^2));
            yy2(i,j) = M1;
        end
 end
 
 % Salva a matriz de parâmetros em um .mat 
[~, filename, ext] = fileparts(file);
savefile = sprintf('%s%s.mat', 'C:\Users\joaop\Google Drive\UFSM\8º Semestre\Filtro Digitais\Trabalho Final\birds\teste\savesteste\', filename);   
save(savefile,'yy1','yy2');

 % Classificação de espécies
 disp ('Iniciando classificação:')
 fprintf('\n');
 
 % Número de arquivos do database
 numfiles = 30;
 
 Corr1 = zeros(1,numfiles);
 
 k = 1;
    for j = 1:numfiles
        clc;
        
        fprintf('Processando algoritmo: %d %% \n',round((j/numfiles)*100));
        % Carrega o arquivo a ser testado
        load(savefile);
        % Parâmetros a ser testado em outras variáveis temporarias
        src1 = yy1;
        src2 = yy2;
        
        % Carregando referência
        disp(sprintf('Carregando database %d... \n',j));
        myfilename2 = sprintf('%sDB%d.mat', 'C:\Users\joaop\Google Drive\UFSM\8º Semestre\Filtro Digitais\Trabalho Final\birds\teste\sample_treino\', j);   ;
        load(myfilename2);
        
        % Parâmetros de referência colocados em variavel temporária
        ref1 = yy1;
        ref2 = yy2;
        
        % Número de linhas de cada matriz
        [timeSrc1, ~] = size(src1);
        [timeRef1, ~] = size(ref1);
        
        % Ajustando o sinal para que o que tenha mais samples seja do mesmo
        % tamanho -> maior tamanho = timeSrc1
        if (timeSrc1<timeRef1)
            src11 = src1;
            src22 = src2;
            clear src1;
            clear src2;
            
            src1 = ref1;
            src2 = ref2;
            clear ref1;
            clear ref2;
            
            ref1 = src11;
            ref2 = src22;
        end
        % Valores ajustados
        [timeSrc, ~] = size(src1);
        [timeRef, ~] = size(ref1);
        
        % Diferença de amostra entre referência e testado
        l = timeSrc-timeRef;
        if(l==0)
            % Se os dois sinais sao do mesmo tamanho, então não é
            % necessário fazer o janelamento
            % Coeficiente de correlação
            corre(1,:) = corr2(src1,ref1); 
            corre(2,:) = corr2(src2,ref2);
            
            % Correlação com cada referência j
            Corr1(k,j) = max(corre); % Pela correlação máxima 
        else
            % Janelamento do sinal a ser testado com a referência
            for i = 1:l
                corre(1,i) = corr2(src1((i:timeRef+i-1),:),ref1);
                corre(2,i) = corr2(src2((i:timeRef+i-1),:),ref2);
            end
            % Correlação máxima 
            maximumCorrelation = zeros(1,2);
           
            for i = 1:2
                
                maximumCorrelation(i) = max(corre(i,:));
                
            end
            Corr1(k,j) = max(maximumCorrelation);     % Máxima correlação com todas as espécies            
        end
        clearvars -except numfiles j k Corr1 filename ext savefile;
    end
     
 
 % Áudios do database (mudar para espécies)
especies =["João-de-Barro","Sabiá-do-Campo","Canário"]; 
nome_imagens = ["joaodebarro","sabiadocampo","canario"];
% Maior valor de correlação e seu respectivo índice
[valor indice] = max(Corr1);

% Threshold para não encontrado no database (arbitrário)
thresh = 7;

% Classifica as espécies
if sum(Corr1) > thresh
if indice>=1 && indice <=10 
    
msg = sprintf('Espécie encontrada! Esse é o canto do(a) %s.',especies(1));
disp(msg);
fileimage = sprintf('%s%s.jpg', 'C:\Users\joaop\Google Drive\UFSM\8º Semestre\Filtro Digitais\Trabalho Final\birds\teste\sample_imagens\', nome_imagens(1));   
figure;
imshow(fileimage);
title (sprintf('%s',especies(1)))

elseif indice >10 && indice <=20
    msg = sprintf('Espécie encontrada! Esse é o canto do(a) %s.',especies(2));
disp(msg);
fileimage = sprintf('%s%s.jpg', 'C:\Users\joaop\Google Drive\UFSM\8º Semestre\Filtro Digitais\Trabalho Final\birds\teste\sample_imagens\', nome_imagens(2));   
figure;
imshow(fileimage);
title (sprintf('%s',especies(2)))

 elseif indice >20 && indice <=30
    msg = sprintf('Espécie encontrada! Esse é o canto do(a) %s.',especies(3));
 disp(msg);
 fileimage = sprintf('%s%s.jpg', 'C:\Users\joaop\Google Drive\UFSM\8º Semestre\Filtro Digitais\Trabalho Final\birds\teste\sample_imagens\', nome_imagens(3));   
figure;
imshow(fileimage);
title (sprintf('%s',especies(3)))
end
else 
    msg = sprintf('Não foi possível identificar o áudio.');
disp(msg);
end
        
        

 
 
 
 