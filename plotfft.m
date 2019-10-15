function [f,mag] = plotfft(y,Fs)
% Plot da FFT de um sinal.
%
%
% Entrada: 
%       y: sinal a ser computado a FFT
%       Fs: Frequência de amostragem.
%       Cor: 'b','r'... cor pra plot do gráfigo
% Saída:
%       f: Frequência 
%       mag: Magnitude do sinal
%              
%
% Autor: João Brum

% FFT do sinal
Y = fft(y);

% Vetor de frequências
w = linspace(0,2*pi,length(Y));
f = w(1:length(w)/2)*(Fs/(2*pi));
mag = abs(Y(1:length(w)/2));

%figure(1)

% Plot do espectro (só metade)
plot(f,mag/length(mag));
xlabel 'Frequência (Hz) '
ylabel 'Módulo'  
grid on