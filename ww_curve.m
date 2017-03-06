% plot Winterwerp curve
figure(1);clf

G = logspace( log10(25),log10(110),20)

KaKb = 3.5e8
C = [0.1; 2; 10; 20]
for i=1:length(C)
   De = KaKb*C(i)./(sqrt(G))
   plot(G,De/1.e6)
   hold on
end
xlim([0 120])
ylim([0 1500])
