%CST- Moisés

clear

%ESTRUTURA
nelem = 160;       %número de elementos da estrutura
nnods = 105;       %número de nós da estrutura
%ELEMENTO CST
dfreedom = 2;    %número de graus de liberdade por nó
noselem = 3;     %número de nós por elemento
% IDENTIFICAÇÃO DOS EIXOS DAS COORDENADAS
x = 0   %eixo e coordenadas na direção x representado por 0
y = 1   %eixo e coordenadas na direção y representado por 1

%PROPRIEDADES DO MATERIAL DE CADA ELEMENTO
for i=1:nelem
E(i) = 200E5    %módulo de elasticidade do material em N/mm^2
v(i) = 0.3     %coeficiente de Poisson
t(i) = 1000    %espessura do cada elemento em mm
end

%GEOMETRIA
%coordenadas do nós
coordenadas =  [0 0; 0 0.5; 0 1; 0 1.5; 0 2; 0.5 0; 0.5 0.5; 0.5 1; 0.5 1.5; 0.5 2; 1 0; 1 0.5; 1 1; 1 1.5; 1 2;1.5 0; 1.5 0.5; 1.5 1; 1.5 1.5; 1.5 2; 2 0; 2 0.5; 2 1; 2 1.5; 2 2; 2.5 0; 2.5 0.5; 2.5 1; 2.5 1.5; 2.5 2; 3 0; 3 0.5; 3 1;3 1.5; 3 2;3.5 0;3.5 0.5;3.5 1;3.5 1.5; 3.5 2; 4 0; 4 0.5; 4 1; 4 1.5; 4 2; 4.5 0; 4.5 0.5; 4.5 1; 4.5 1.5; 4.5 2; 5 0; 5 0.5; 5 1; 5 1.5; 5 2; 5.5 0; 5.5 0.5; 5.5 1; 5.5 1.5; 5.5 2; 6 0; 6 0.5; 6 1; 6 1.5; 6 2; 6.5 0; 6.5 0.5; 6.5 1; 6.5 1.5; 6.5 2; 7 0; 7 0.5; 7 1; 7 1.5; 7 2; 7.5 0; 7.5 0.5; 7.5 1; 7.5 1.5; 7.5 2; 8 0; 8 0.5; 8 1; 8 1.5; 8 2; 8.5 0; 8.5 0.5; 8.5 1; 8.5 1.5; 8.5 2; 9 0; 9 0.5; 9 1; 9 1.5; 9 2; 9.5 0; 9.5 0.5; 9.5 1; 9.5 1.5; 9.5 2; 10 0; 10 0.5; 10 1; 10 1.5; 10 2];   
%conectividades do elemento
conectividades =  [1 6 7; 1 7 2; 2 7 8; 2 8 3; 3 8 4; 8 9 4; 4 9 5; 9 10 5; 6 11 12; 6 12 7; 7 12 13; 7 13 8; 8 13 9; 13 14 9; 9 14 10; 14 15 10; 11 16 17; 11 17 12; 12 17 18; 12 18 13; 13 18 14; 18 19 14; 14 19 15; 19 20 15; 16 21 22; 16 22 17; 17 22 23; 17 23 18; 18 23 19; 23 24 19; 19 24 20; 24 25 20; 21 26 27; 21 27 22; 22 27 28; 22 28 23; 23 28 24; 28 29 24; 24 29 25; 29 30 25; 26 31 32; 26 32 27; 27 32 33; 27 33 28; 28 33 29; 33 34 29; 29 34 30; 34 35 30; 31 36 37; 31 37 32; 32 37 38; 32 38 33; 33 38 34; 38 39 34; 34 39 35; 39 40 35; 36 41 42; 36 42 37; 37 42 43; 37 43 38; 38 43 39; 43 44 39; 39 44 40; 44 45 40; 41 46 47; 41 47 42; 42 47 48; 42 48 43; 43 48 44; 48 49 44; 44 49 45; 49 50 45; 46 51 52; 46 52 47; 47 52 53; 47 53 48; 48 53 49; 53 54 49; 49 54 50; 54 55 50; 51 56 57; 51 57 52; 52 57 58; 52 58 53; 53 58 54; 58 59 54; 54 59 55; 59 60 55; 56 61 62; 56 62 57; 57 62 63; 57 63 58; 58 63 59; 63 64 59; 59 64 60; 64 65 60; 61 66 67; 61 67 62; 62 67 68; 62 68 63;63 68 64; 68 69 64; 64 69 65; 69 70 65; 66 71 72; 66 72 67; 67 72 73; 67 73 68; 68 73 69; 73 74 69; 69 74 70; 74 75 70; 71 76 77; 71 77 72; 72 77 78; 72 78 73; 73 78 74; 78 79 74; 74 79 75; 79 80 75; 76 81 82; 76 82 77; 77 82 83; 77 83 78; 78 83 79; 83 84 79; 79 84 80; 84 85 80; 81 86 87; 81 87 82; 82 87 88; 82 88 83; 83 88 84; 88 89 84; 84 89 85; 89 90 85; 86 91 92; 86 92 87; 87 92 93; 87 93 88; 88 93 89; 93 94 89; 89 94 90; 94 95 90; 91 96 97; 91 97 92; 92 97 98; 92 98 93; 93 98 94; 98 99 94; 94 99 95; 99 100 95; 96 101 102; 96 102 97; 97 102 103; 97 103 98; 98 103 99; 103 104 99; 99 104 100; 104 105 100]

    %CONDIÇÕES DE CONTORNO

ncargas = 1;                                    %número de cargas aplicadas na estrutura
nrestricao = 10;                               %número de restrições na estrutura
carga = [103 y -10000];                       %carga concentrada
restricao = [1 x; 1 y; 2 x; 2 y; 3 x; 3 y; 4 x; 4 y; 5 x; 5 y];   %localização e direção da restrição

%CÁLCULO DA MATRIZ CONSTITUTIVA DE CADA ELEMENTO
for k = 1:nelem
    C =(E(k)/(1-v(k)^2))*[1 v(k) 0; v(k) 1 0; 0 0 ((1-v(k))/2)];    %matriz constitutiva

%MATRIZ RIGIDEZ E SUPERPOSIÇÃO

for i = 1:nnods*dfreedom
        Lglobal(i) = 0;   %vetor de cargas
    for j = 1:nnods*dfreedom
        Kglobal(i,j) = 0;  %matriz rigidez global
    end
end

for i = 1:nelem
    %coordenadas dos pontos do elemento
    for j = 1:noselem
        X(i,j) = coordenadas(conectividades(i,j),1);
        Y(i,j) = coordenadas(conectividades(i,j),2);
    end
%CÁLCULO DA ÁREA DO ELEMENTO          
    A(i) = 0.5*det([1 X(i,1) Y(i,1); 1 X(i,2) Y(i,2); 1 X(i,3) Y(i,3)]);

%CÁLCULO DOS COEFICIENTES DE FORMA

a(i,1) = X(i,2)*Y(i,3) - X(i,3)*Y(i,2);
a(i,2) = X(i,3)*Y(i,1) - X(i,1)*Y(i,3);
a(i,3) = X(i,1)*Y(i,2) - X(i,2)*Y(i,1);
b(i,1) = Y(i,2) - Y(i,3);
b(i,2) = Y(i,3) - Y(i,1);
b(i,3) = Y(i,1) - Y(i,2);
c(i,1) = X(i,3) - X(i,2);
c(i,2) = X(i,1) - X(i,3);
c(i,3) = X(i,2) - X(i,1);
end

for i = 1:nelem

%CÁLCULO DA MATRIZ DE APROXIMAÇÃO DE DEFORMAÇÃO

B = (1/(2*A(i)))*[b(i,1) 0 b(i,2) 0 b(i,3) 0;0 c(i,1) 0 c(i,2) 0 c(i,3); c(i,1) b(i,1) c(i,2) b(i,2) c(i,3) b(i,3)];
%cálcula-se através da multiplicação entre as matrizes B e C

%MATRIZ RIGIDEZ

K = A(i)*t(i)*B'*C*B;

%SUPERPOSIÇÃO

for j = 1:noselem
    for k = 1:dfreedom
        colE = (j-1)*dfreedom + k;
        for l = 1:noselem
            for m = 1:dfreedom
                linE = (l-1)*dfreedom + m;
                colG = (conectividades(i,j)-1)*dfreedom + k;
                linG = (conectividades(i,l)-1)*dfreedom + m;
                Kglobal(linG,colG) = Kglobal(linG,colG) + K(linE,colE);
            end
        end
    end
end
end
end

%MATRIZ CARREGAMENTO

for i = 1:ncargas
    linG = 2*carga(i,1)-(1-carga(i,2));
    Lglobal(linG) = Lglobal(linG) + carga(i,3);
end

%CONDIÇÕES DE CONTORNO

for i = 1:nrestricao
    linG = 2*(restricao(i,1))-(1-(restricao(i,2)));
    Lglobal(linG)=0
    for j = 1:nnods*dfreedom
        if linG == j;
            Kglobal(linG,j) = 1;
        else Kglobal(linG,j) = 0;
        end
        if linG == j;
            Kglobal(j,linG) = 1;
        else Kglobal(j,linG) = 0;
        end
    end
end 

%DESLOCAMENTOS

desloc = Lglobal/Kglobal

%Tensoes em cada elemento

for i = 1:nelem
    C = (E(i)/(1-v(i)^2))*[1 v(i) 0; v(i) 1 0; 0 0 ((1-v(i))/2)];
    B = (1/(2*A(i)))*[b(i,1) 0 b(i,2) 0 b(i,3) 0;0 c(i,1) 0 c(i,2) 0 c(i,3); c(i,1) b(i,1) c(i,2) b(i,2) c(i,3) b(i,3)];
    desloc1 = [desloc(2*conectividades(i,1) - 1); desloc(2*conectividades(i,1)); desloc(2*conectividades(i,2) - 1); desloc(2*conectividades(i,2)); desloc(2*conectividades(i,3) - 1); desloc(2*conectividades(i,3))];
    Tensao = C*B*desloc1
   
end

%Fator para visualização da deformada
distmax(1,1)=max(coordenadas(:,1)) - min(coordenadas(:,1));
distmax(2,1)=max(coordenadas(:,2)) - min(coordenadas(:,2));

deslabsoluto=abs(desloc)

f=(0.15*max(distmax(:,1)))/(max(deslabsoluto))

desloc
axis off 
daspect([1,1,1])
hold on;

for i = 1:nelem
   no1 = conectividades(i,1);
   no2 = conectividades(i,2);
   no3 = conectividades(i,3);
   x(1) = coordenadas(no1,1); y(1) = coordenadas(no1,2);
   x(2) = coordenadas(no2,1); y(2) = coordenadas(no2,2);
   x(3) = coordenadas(no3,1); y(3) = coordenadas(no3,2);
   x(4) = x(1);               y(4) = y(1);    
   plot(x,y,'LineWidth',3)
   plot(x,y,'o')
end
%coordenadas dos nós deformados
for i = 1:nnods
        defx(i) = coordenadas(i,1)+ f*desloc(2*i-1);
        defy(i) = coordenadas(i,2)+ f*desloc(2*i);
end

for i = 1:nnods
    deformacao(i,1) = defx(i);
    deformacao(i,2) = defy(i);
end

for i = 1:nelem
   no1 = conectividades(i,1);
   no2 = conectividades(i,2);
   no3 = conectividades(i,3);
   xdef(1) = deformacao(no1,1); ydef(1) = deformacao(no1,2);
   xdef(2) = deformacao(no2,1); ydef(2) = deformacao(no2,2);
   xdef(3) = deformacao(no3,1); ydef(3) = deformacao(no3,2);

   plot(xdef,ydef,'r','LineWidth',1)
   plot(xdef,ydef,'o')
end

hold off;   