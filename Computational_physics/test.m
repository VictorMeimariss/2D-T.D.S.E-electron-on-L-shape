clc;
clear;

% Συντελεστές αριθμητή και παρονομαστή της G(s)F(s)
K = 1;  % Δεν επηρεάζει το root locus, μπορείς να το αφήσεις ως 1

% G(s)F(s) = K*(s+5) / (s*(s+10)*(s+3))
numerator = K * [1 5];                   % (s + 5)
denominator = conv([1 0], conv([1 10], [1 3]));  % s*(s+10)*(s+3)

% Δημιουργία μεταβλητής μεταφοράς
sys = tf(numerator, denominator);

% Σχεδίαση γεωμετρικού τόπου των ριζών
figure;
rlocus(sys)
title('Γεωμετρικός Τόπος των Ριζών (Root Locus)')
grid on;