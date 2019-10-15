function [signals1,n] = filterbank(y,Fs,t,n_filtros)
% Design de tamanho de n_filtros filtros passa banda, dividindo o spectro em
% sub bandas + filtragem de sinal (sinal filtrado pelas sub-bandas).
%
%
% Entrada: 
%       y: sinal a ser filtrado
%       Fs: Frequência de amostragem
%       t: vetor de tempo
%       n_filtros: número de filtros
%
% Saída:
%       signals1: sinal filtrado (cada coluna é uma banda)
%       n: numero de filtros
%      
%
%Autor: João Brum
%

% Comprimento do sinal a ser filtrado
l = length(y);
% Calcula a FFT
A = abs(fft(y));
% Vetor de frequência correspondente
f = linspace(0,2*pi,length(A))*(Fs/(2*pi));

% Banda total Fs/2 - 100 (começa em 100)
bandatotal = floor(f(end)/2) - 100;
% Tamanho da banda
banda = floor(bandatotal/n_filtros);
% Armazena o sinal filtrado por cada banda
signals1 = zeros(l,n_filtros);
%figure(3); subplot(211);hold on;
for n = 1:n_filtros
 % Frequência de início e fim 
start1 = ((n-1)*banda)+100;
stop1 = ((n)*banda);

% Filtro passa banda butterworth, de ordem 5 
[b,a] =butter(5,[start1/(Fs/2) stop1/(Fs/2)],'bandpass');
% Filtra o sinal
yfilt = filter(b,a,y);

% Armazena o sinal 
signals1(:,n)= yfilt;
% Resposta em frequência do filtro
[H,f] = freqz(b,a,1024,Fs);

%Plots
%plot(f,20*log10(abs(H)));
%plot(t,yfilt);
%plotfft(yfilt,Fs);
legendInfo{n} = ['Banda Nº - ' num2str(n)]; 

end
% title('Sinal filtrado para cada banda');
% xlabel('Tempo (s)');
% ylabel('Amplitude (V)');
% legend(legendInfo,'Location','eastoutside','NumColumns',2);
% grid on

% subplot(212);
% plot(t,y);
% title('Sinal Original');
% xlabel('Tempo (s)');
% ylabel('Amplitude (V)');
% grid on

end

