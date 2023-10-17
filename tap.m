clear all;
x_init = linspace(-10,10,1001);
y_init = sinc(1.5*x_init);
figure(1);
subplot(2,1,1);
plot(x_init,y_init);
grid on;
x_sample = linspace(-10, 10, 21);
y_sample = [];
for i=1:50:1001
y_sample = [y_sample y_init(i)];
end
figure(1);
subplot(2,1,2);
stem(x_sample,y_sample);
grid on;
y_input_row = [0 0 0 0 0 0 0 0 0 0 fliplr(y_sample) 0 0 0 0 0 0 0 0 0 0];
Pmatrix = [];
Peq_matrix = [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
for i=21:41
rowx = fliplr(y_input_row(i-20:i));
Pmatrix = [Pmatrix; rowx];
end
Cmatrix = inv(Pmatrix)*Peq_matrix';
Pfinal = Pmatrix*Cmatrix;
y_final = sinc(x_init);
j = 1;
for i=1:50:1001
y_final(i) = y_init(i)*Pfinal(j);
j = j+1;
end
figure(2);
plot(x_init,y_final);
hold on;
stem(x_sample,Pfinal');
grid on;