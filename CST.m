%CST- Moisés

clear

%ESTRUTURA
nelem = 2;       %número de elementos da estrutura
nnods = 4;       %número de nós da estrutura
%ELEMENTO CST
dfreedom = 2;    %número de graus de liberdade por nó
noselem = 3;     %número de nós por elemento
% IDENTIFICAÇÃO DOS EIXOS DAS COORDENADAS
x = 0   %eixo e coordenadas na direção x representado por 0
y = 1   %eixo e coordenadas na direção y representado por 1

%PROPRIEDADES DO MATERIAL DE CADA ELEMENTO
for i=1:nelem
E(i) = 1    %módulo de elasticidade do material em N/mm^2
v(i) = 0.000001     %coeficiente de Poisson
t(i) = 1    %espessura do cada elemento em mm
end

rho=1;

thic=1;



%GEOMETRIA
%coordenadas do nós
coordenadas =  [0 0; 1 0; 1 1; 0 1];   
%conectividades do elemento
conectividades =  [1 2 3; 1 3 4]

    %CONDIÇÕES DE CONTORNO

ncargas = 2;                                    %número de cargas aplicadas na estrutura
nrestricao = 4;                               %número de restrições na estrutura
carga = [2 x 0.5;3 x 0.5];                       %carga concentrada
restricao = [1 x; 1 y; 4 x; 4 y];   %localização e direção da restrição

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
    
    
    M=rho*A(i)*thic/3*eye(3);

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