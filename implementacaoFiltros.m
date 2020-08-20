#[amostras, frequência de amostragem] 400312, 40000
[a,fs] = audioread('C:\Users\meliodas\Desktop\5Semestre\PDS\work\xan.wav'); #Audio representando xa[n]
[b,fs] = audioread('C:\Users\meliodas\Desktop\5Semestre\PDS\work\xbn.wav'); #Audio representando xb[n]

#linspace(início, fim, n° de amostras em a)
ta = linspace(0,length(a)/fs,length(a)); #Vetor do tempo total em segundos
tb = linspace(0,length(b)/fs,length(b)); #Vetor do tempo total em segundos

Nfft = 65536; #Tamanho da fft

fv = linspace(-20000,20000,Nfft); #Vetor da Frequência Normalizado

A=fftshift(abs(fft(a,Nfft))); #Vetor de frequência
B=fftshift(abs(fft(b,Nfft))); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DOWNSAMPLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####FILTRO IIR Butterworth####
wp = 4000/(fs/2);
ws = 8000/(fs/2);
[n,wn]=buttord(wp,ws,3,35);
[p,o]=butter(n,wn); #Retorna ordem e frequência de corte normalizada
fa=filter(p,o,a); #Sinal xa[n] filtrado
fb=filter(p,o,b); #Sinal xb[n] filtrado
FA = fftshift(abs(fft(fa,Nfft))); #Espectro de amplitude de xa[n] filtrado
FB = fftshift(abs(fft(fb,Nfft))); #Espectro de amplitude de xa[n] filtrado
figure(1)
freqz(p,o) #Função de Transferência do Filtro

figure(2)
subplot(221)
plot(ta,a)
title('xa[n]')
subplot(222)
plot(fv,A)
title('Espectro Amplitude xa[n]')
subplot(223)
plot(ta,fa)
title('xa[n] filtrado')
subplot(224)
plot(fv, FA)
title('Espectro Amplitude xa[n] filtrado')
#sound(fa,fs)

figure(3)
subplot(221)
plot(tb,b)
title('xb[n]')
subplot(222)
plot(fv,B)
title('Espectro Amplitude xb[n]')
subplot(223)
plot(tb,fb)
title('xb[n] filtrado')
subplot(224)
plot(fv, FB)
title('Espectro Amplitude xb[n] filtrado')
#sound(fb,fs)

####Dizimação####
M = 5; #Fator de Dizimação
##Para xa[n]
for i=1:length(a)/M
  aD(i) = fa(i*5);   #y[n] = x[nM]
end
for i=1:length(b)/M
  bD(i) = fb(i*5);
end  
fsd = length(aD)/(length(a)/fs); #8Hz obtido após downsampling, é igual para bD
NfftD = 8192; #Novo tamanho da fft de acordo com a nova frequência de amostragem obtida
tda=linspace(0,length(aD)/fsd,length(aD)); #Vetor de tempo total em segundos
tdb=linspace(0,length(bD)/fsd,length(bD)); #Vetor de tempo total em segundos
fvD = linspace(-4000,4000,NfftD); #Vetor da frequência normalizado com nova Nfft
Ad=fftshift(abs(fft(aD,NfftD))); #Espectro de amplitude do xa[n] dizimado
Bd=fftshift(abs(fft(bD,NfftD))); #Espectro de amplitude do xa[n] dizimado
figure(4)
subplot(221)
plot(ta,a)
title('xa[n]')
xlabel('Tempo')
ylabel('Amplitude')
subplot(222)
plot(fv,A)
title('Espectro Amplitude xa[n]')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
subplot(223)
plot(tda,aD)
title('xa[n] Downsampling')
xlabel('Tempo')
ylabel('Amplitude')
subplot(224)
plot(fvD, Ad)
title('Espectro Amplitude xa[n] Downsampling')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
#sound(aD,fsd)
#stem(aD);

figure(5)
subplot(221)
plot(tb,b)
title('xb[n]')
xlabel('Tempo')
ylabel('Amplitude')
subplot(222)
plot(fv,B)
title('Espectro Amplitude xb[n]')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
subplot(223)
plot(tdb,bD)
title('xb[n] Downsampling')
xlabel('Tempo')
ylabel('Amplitude')
subplot(224)
plot(fvD, Bd)
title('Espectro Amplitude xb[n] Downsampling')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
#sound(bD,fsd)
#stem(bD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODULAÇÃO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

####Modulação AM DSB####
t1=0:1/fsd:length(tda)/fsd - 1/fsd;
t2=0:1/fsd:length(tda)/fsd - 1/fsd;
##Para xa[n]
f1=2E3;
c1=2*cos(2*pi*f1*t1);
aDM=(aD).*c1;
ADM=fftshift(abs(fft(aDM,NfftD)));
figure(6)
subplot(221)
plot(tda,aD)
title('xaD[n] Downsampling')
xlabel('Tempo')
ylabel('Amplitude')
subplot(222)
plot(fvD,Ad)
title('Espectro Amplitude xaD[n] Downsampling')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
subplot(223)
plot(t1,aDM)
title('xaD[n] Downsampling Modulado')
xlabel('Tempo')
ylabel('Amplitude')
subplot(224)
plot(fvD, ADM)
title('Espectro Amplitude xaD[n] Downsampling Modulado')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
#sound(aDM,fsd)

##Para xb[n]
f2=1E3;
c2=cos(2*pi*f2*t1);
bDM=(bD).*c2;
BDM=fftshift(abs(fft(bDM,NfftD)));
figure(7)
subplot(221)
plot(tdb,bD)
title('xbD[n] Downsampling')
xlabel('Tempo')
ylabel('Amplitude')
subplot(222)
plot(fvD,Bd)
title('Espectro Amplitude xbD[n] Downsampling')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
subplot(223)
plot(t1,bDM)
title('xbD[n] Downsampling Modulado')
xlabel('Tempo')
ylabel('Amplitude')
subplot(224)
plot(fvD, BDM)
title('Espectro Amplitude xbD[n] Downsampling Modulado')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
#sound(bDM,fsd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOMA DOS SINAIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss = aDM + bDM;
SS=fftshift(abs(fft(ss,NfftD)));
figure(8)
subplot(211)
plot(tda, ss)
title('aDM + bDM')
xlabel('Tempo')
ylabel('Amplitude')
subplot(212)
plot(fvD, SS)
title('Espectro de amplitude de aDM + bDM')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
sound(ss,fsd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTRO PASSA FAIXA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = 3;
ks = 60;
fs1 = 1000;
fp1 = 1500;
fp2 = 3000;
fs2 = 3500;

#Frequências normalizadas
wp1 = fp1/(fsd/2);
wp2 = fp2/(fsd/2);
ws1 = fs1/(fsd/2);
ws2 = fs2/(fsd/2);
[N wc] = buttord([wp1 wp2], [ws1 ws2], kp, ks);
[h j] = butter(N, wc, 'bandpass'); 
pa = filter(h,j,ss); #Isola canal 1, xa[n]
PA = fftshift(abs(fft(pa,NfftD)));
figure(9)
freqz(b,a);
figure(10)
subplot(211)
plot(tda, pa)
title('Após Filtro Passa Faixa')
xlabel('Tempo')
ylabel('Amplitude')
subplot(212)
plot(fvD, PA)
title('Espectro: Após Filtro Passa Faixa')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
#sound(pa,fsd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEMODULAÇÃO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dm = pa.*c1;
DM = fftshift(abs(fft(dm,NfftD)));
figure(11)
subplot(211)
plot(tda,dm)
title('Sinal Demodulado')
xlabel('Tempo')
ylabel('Amplitude')
subplot(212)
plot(fvD,DM)
title('Espectro:Sinal Demodulado')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
#sound(dm,fsd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPSAMPLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####EXPANSÃO####
L = M; #Fator de expansão L = M = 5;
for s=1:length(a)/L
  adE(s)=0;
end

for s=1:length(a)/L
  adE(s*5)=dm(s);
end
tadE=linspace(0,length(adE)/fs,length(adE));
ADE = fftshift(abs(fft(adE,Nfft)));
figure(12)
subplot(211)
plot(tadE, adE)
title('adE[n]')
xlabel('Tempo')
ylabel('Amplitude')
subplot(212)
plot(fv, ADE)
title('Espectro de magnitude adE[n]')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
#stem(adE)

####FILTRO IIR Butterworth FINAL####
wpp = 1200/(fs/2);
wss = 4530/(fs/2);
[nn,wnn]=buttord(wpp,wss,3,80);
[pp,oo]=butter(nn,wnn);
y=filter(pp,oo,adE); #Sinal filtrado, tentativa de reconstrução de xa[n]
Y = fftshift(abs(fft(y,Nfft))); #Espectro de amplitude de sinal filtrado
y =(max(a)/max(y))*y;
Y =(max(A)/max(Y))*Y;
ty=linspace(0,length(y)/fs,length(y)); #Vetor do tempo total em segundos
figure(13)
subplot(221)
plot(ty,y)
title('y[n]')
xlabel('Tempo')
ylabel('Amplitude')
subplot(223)
plot(fv, Y)
title('Espectro de amplitude y[n]')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
subplot(222)
plot(ta,a)
title('xa[n]')
xlabel('Tempo')
ylabel('Amplitude')
subplot(224)
plot(fv,A)
title('Espectro de amplitude xa[n]')
xlabel('Frequência(Hz)')
ylabel('Amplitude')
#sound(a,fs);
#sound(y,fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%& Análise de PSNR e Erro Quadrático &&&&&&&&&&&&&&&&&&&&&&&&&
y(400311)=0;
y(400312)=0;
err = immse(a',y);
fprintf('MSE é :%0.5f\n',err);
[peaksnr snr] = psnr(y,a');
fprintf('O PSNR é :%0.5f\n',peaksnr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%