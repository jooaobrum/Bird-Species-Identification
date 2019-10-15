% Reconhecimento de espécies de pássaros através do canto: Estágio de
% Treinamento.
%
% Como funciona o database:
% DB1.mat-DB10.mat = João-de-Barro
% DB11.mat-DB20.mat = Sabiá-do-Campo
% DB21.mat-DB30.mat = Canário
% 
% Nota: Manter o database atualizado (numfiles LINHA 19)
% Autor: João Brum

clc
clear all
close all


warning('off');

numfiles = 30;

for k = 1:numfiles
% Seleciona o arquivo na pasta para testar o canto
file = sprintf('%sDB%d.wav', 'C:\Users\joaop\Google Drive\UFSM\8º Semestre\Filtro Digitais\Trabalho Final\birds\teste\sample_treino\', k);   
 
% Armazena o áudio para treinamento (5 segundos equivalente em amostras)
[y,Fs] = audioread(file);

% Divide o áudio em parte com "voz" e sem "voz", agrupando todas as com
% "voz"
segmentos = detectVoiced(y,Fs);
clear y;
y = segmentos;

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

%% Plot do dominio tempo depois da decimação
t = 0:1/Fsnew:(length(Y1)-1)*(1/Fsnew);
plot(t,Y1);
xlabel 'Tempo(s)'
ylabel 'Amplitude (V)'

%% Plot do dominio frequência 
plotfft(Y1,Fsnew);

%% Divide o sinal no domínio tempo em 40ms
frame = round(0.04*Fsnew);

% Número de segmentos
segmentos = floor(length(Y1)/frame);

for i = 1:segmentos
    Y2 = Y1((i-1)*frame+1:i*frame);
    F1 = abs(fft(Y2));
    F1 = F1(1:length(F1)/2);
    f = linspace(0,pi*(Fsnew/(2*pi)),length(F1));
    [val ind] = max(F1);
    Freq1 = f(ind);
end


%% Divide o sinal para passar pelo banco de filtros

% Número de bandas 
n_filtros = 30;
[freqSignals,n] = filterbank(Y1,Fsnew,t,n_filtros);

% Armazena as matrizes de Parâmetros
% yy1 -> componente da frequência dominante de cada segmento de cada banda
yy1 = zeros(segmentos,n);
% yy2 -> energia média de cada segmento de cada banda
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
 
 % Salva os arquivos em .mat
 savefile = sprintf('%sDB%d.mat', 'C:\Users\joaop\Google Drive\UFSM\8º Semestre\Filtro Digitais\Trabalho Final\birds\teste\sample_treino\', k);   
  
 save(savefile,'yy1','yy2');
 
 clearvars -except k
end

