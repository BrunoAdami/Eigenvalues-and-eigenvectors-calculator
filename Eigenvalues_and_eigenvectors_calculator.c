//
//  main.c
//  EP1_calculo_numerico
//
//  Created by Bruno Adami Serine on 1/25/19.
//  Copyright © 2019 Bruno Adami. All rights reserved.
//

#include <stdio.h>
#include <math.h> /*essa biblioteca foi incluida para ser utilizada a funcao sqrt() de raiz quadrada, 
uma operação que já é nativa em Python, mas não em C */

#define MAX 100 //utilizado como tamanho máximo para as matrizes e vetores

//a funcao imprimeMatriz imprime uma matriz quadrada

void imprimeMatriz (double matriz[MAX][MAX], int tamanho){
    for (int i = 0; i < tamanho; i++){
        for (int j = 0; j < tamanho; j++){
            printf("%lf ", matriz[i][j]);
            
            if (j + 1 == tamanho) printf("\n");
        }
    }
    printf("\n");
}

// a funcao imprimeAutovaloresAutovetores recebe as matrizes de autovalores,
// onde os autovalores estão em sua diagonal principal, e de autovetores, onde
// os autovetores estão nas colunas desta matriz, e imprime o autovalor seguido
// de seu autovetor associado

void imprimeAutovaloresAutovetores (double autovalores[MAX][MAX], double autovetores[MAX][MAX],
                                    int tamanho) {
    printf("\nAssociated eigenvalues ​​and eigenvectors:\n\n");
    
    for (int a = 0; a < tamanho; a++){
        printf("Eigenvalue %d = %lf  |  ", a + 1, autovalores[a][a]);
        printf("Eigenvector %d = (", a + 1);
        
            for (int b = 0; b < tamanho; b++){
                if (b == tamanho - 1) printf("%lf);\n", autovetores[b][a]);
                else printf("%lf, ", autovetores[b][a]);
        }
    }
}

// a funcao imprimeVetorVertical imprime o vetor armazenado na primeira coluna da matriz vetor

void imprimeVetorVertical (double vetor[MAX][MAX], int tamanhoDoVetor){
    for (int i = 0; i < tamanhoDoVetor; i++)
        printf("%lf\n", vetor[i][0]);
    printf("\n");
}

// funcao calculaNorma calcula a norma do vetor escrito na coluna indicada da matriz de vetores

double calculaNorma (double vetores[MAX][MAX], int tamanho, int coluna){
    double norma = 0;

    for (int i = 0; i < tamanho; i++){
        norma += (vetores[i][coluna - 1] * vetores[i][coluna - 1]); //soma as coordenadas do vetor ao quadrado
    }

    norma = sqrt(norma); 
    return norma;

}

// a funcao geraIdentidade produz uma matriz identidade de tamanho n

void geraIdentidade (double identidade[MAX][MAX], int tamanho){
    for (int i = 0; i < tamanho; i++){
        for (int j = 0; j < tamanho; j++){
            if (i == j) identidade [i][j] = 1;
            else identidade[i][j] = 0;
        }
    }
}

// funcao verificaSinal retorna -1 se o número passado for negativo, +1 caso contrário

int verificaSinal (double valor){
    if (valor < 0) return -1;
    else return 1;
}

void transpoeMatriz (double original[MAX][MAX], double transposta[MAX][MAX], int tamanho){
    for (int i = 0; i < tamanho; i++){
        for (int j = 0; j < tamanho; j++){
            transposta[i][j] = original[j][i];
        }
    }
}

// funcao somaColunas soma os vetores de duas matrizes presentes em uma determinada coluna,
// colocando o vetor resultado na primeira coluna da matriz soma

void somaColunas (double matriz1[MAX][MAX], double matriz2[MAX][MAX], double soma[MAX][MAX], int tamanho, int coluna){
    for (int i = 0; i < tamanho; i++) 
        soma[i][0] = matriz1[i][coluna - 1] + matriz2[i][coluna - 1];
}

//subtraiMatrizes faz a subtracao matriz1 - matriz2 e armazena em resultado

void subtraiMatrizes (double matriz1[MAX][MAX], double matriz2[MAX][MAX], double resultado[MAX][MAX], int tamanho){
    for (int i = 0; i < tamanho; i++){
        for (int j = 0; j < tamanho; j++) resultado[i][j] = matriz1[i][j] - matriz2[i][j];
    }
}

void multiplicaMatrizPorConstante (double matriz[MAX][MAX], double constante, int tamanho){
    for (int i = 0; i < tamanho; i++){
        for (int j = 0; j < tamanho; j++){
            matriz[i][j] *= constante;
        }
    }
}

// a funcao copiaMatriz faz uma copia da matriz original na matriz copia (as matrizes precisam ser quadradas)

void copiaMatriz (double original[MAX][MAX], double copia[MAX][MAX], int tamanho){
    for (int i = 0; i < tamanho; i++){
        for (int j = 0; j < tamanho; j++) copia[i][j] = original[i][j];
    }
}

// produtoTensorialDeV calcula o produto tensorial do vetor V armazenado na primeira coluna da 
// matriz vetor por ele mesmo, armazenando a matriz gerada em resultado

void produtoTensorialDeV (double vetor[MAX][MAX], double resultado[MAX][MAX], int tamanho){
        for (int a = 0; a < tamanho; a++){
            for (int b = 0; b < tamanho; b++){
                resultado[a][b] = vetor[a][0] * vetor[b][0];
            }
        }
}

//produtoInternoDeV calcula e retorna o produto interno do vetor V armazanado na primeira
// coluna da matriz vetor

double produtoInternoDeV (double vetor[MAX][MAX], int tamanho){
    double produtoInt = 0;

    for (int i = 0; i < tamanho; i++) produtoInt += (vetor[i][0] * vetor[i][0]);

    return produtoInt;
}

// zeraAcimaDaDiagonal zera os elementos da matriz passada que estão acima da diagonal principal

void zeraAcimaDaDiagonal (double matriz[MAX][MAX], int tamanho){
    for (int i = 0; i < tamanho - 1; i++){
        for (int j = i + 1; j < tamanho; j++) matriz[i][j] = 0;
    }
}

//funcao geraVn produz o vetor Vn que será utilizado para gerar a transformação de Householder
//esta funcao nao altera a matriz A passada como parâmetro

void geraVn (double A[MAX][MAX], double Vn[MAX][MAX], int tamanho, int n){

    double identidade[MAX][MAX], A2[MAX][MAX]; //será uma matriz identidade, e portanto terá os vetores da base canônica de R^n em suas colunas
    
    copiaMatriz(A, A2, tamanho);
    zeraAcimaDaDiagonal(A2,tamanho);
    geraIdentidade(identidade, tamanho);
    multiplicaMatrizPorConstante(identidade, ( verificaSinal(A2[n - 1][n - 1]) * calculaNorma(A2, tamanho, n) ), tamanho);
    somaColunas(A2, identidade, Vn, tamanho, n); //armazena o vetor Vn gerado na primeira coluna da matriz Vn
}

 // geraHouseholder produz a matriz da transformação de Householder em H[MAX][MAX]
 // a partir do vetor aramazenado na primeira coluna da matriz vetor[MAX][MAX]

void geraHouseholder( double H[MAX][MAX], double vetor[MAX][MAX], int tamanho){

    double constante = 2 * (1/produtoInternoDeV(vetor, tamanho));
    double identidade[MAX][MAX], produtoExterno[MAX][MAX];

    geraIdentidade(identidade, tamanho);
    produtoTensorialDeV(vetor, produtoExterno, tamanho);
    multiplicaMatrizPorConstante(produtoExterno, constante, tamanho);
    subtraiMatrizes(identidade, produtoExterno, H, tamanho);

}


void transpoeVetor (double vetor[MAX][MAX], double vetorTransposto[MAX][MAX], int tamanhoVetor){

    int n = tamanhoVetor;
    for (int i = 0; i < n; i++){
        vetorTransposto[i][0] = vetor[0][i];
    }
    
}


// a funcao lerMatriz solicita uma matriz ao usuario elemento por elemento e armazena no programa a matriz
// a funcao tb armazena o tamanho da matriz em n

void lerMatriz (double matriz[MAX][MAX], int *n){
    int tamanho, leMatriz = 1, respostaValida ;
    char resposta;
    
    while (leMatriz == 1){
        respostaValida = 0;
        printf("Forneça a ordem da matriz quadrada A:\n");
        scanf("%d", &tamanho);
        printf("Forneça agora os elementos da matriz A.\n");
    
        for (int i = 0; i < tamanho; i++){
            for (int j = 0; j < tamanho; j++){
                printf("a%d%d = ", i+1, j+1 );
                scanf("%lf", &matriz[i][j]);
            }
        }
    
        *n = tamanho;
        printf("\nMatriz armazenada:\n\n");
        imprimeMatriz(matriz, tamanho);
        while (respostaValida == 0){
                printf("Sua matriz está correta? [s/n]: ");
                scanf(" %c", &resposta);
            
                if (resposta == 's') {
                    leMatriz = 0;
                    respostaValida = 1;
                    printf("\n");
                } else if (resposta == 'n') {
                    respostaValida = 1;
                    printf("\nTudo bem, vamos ler sua matriz novamente!\n");
                } else printf("\nResposta inválida. Responda apenas com \"s\" ou \"n\".\n");
        }
    }
    
}

// a funcao produtoMatrizes faz o produto de duas matrizes quadradas de tamanho n

void produtoMatrizes (double matriz1[MAX][MAX], double matriz2[MAX][MAX], double produto[MAX][MAX], int tamanho){
    double soma = 0;

    for (int i = 0; i < tamanho; ++i){
        for (int j = 0; j < tamanho; ++j){
            for (int k = 0; k < tamanho; ++k){
                soma += matriz1[i][k]*matriz2[k][j];
            }

            produto[i][j] = soma;
            soma = 0;
        }
    }
}

//fatoraQR executa a fatoração OR da matriz quadrada A

void fatoraQR (double A[MAX][MAX], double Q[MAX][MAX],double R[MAX][MAX], int tamanho){

    double Hvn[MAX][MAX]; //matriz de householder para o vetor Vn
    double Vn[MAX][MAX]; //matriz que irá armazenar Vn em sua primeira coluna
    double Qparcial[MAX][MAX]; //usado na função do produto
    double Rparcial[MAX][MAX];

    geraIdentidade(Qparcial, tamanho);

    for (int n = 1; n <= tamanho - 1; n++){

        if (n == 1) geraVn(A, Vn, tamanho, n); 
        else geraVn(Rparcial, Vn, tamanho, n);
        
        geraHouseholder(Hvn, Vn, tamanho);
        produtoMatrizes(Qparcial, Hvn, Q, tamanho);
        copiaMatriz(Q, Qparcial, tamanho); //copio o resultado do atual produto para a matriz Qparcial para usar no próximo produto
        
        if (n == 1){
                produtoMatrizes(Hvn, A, R, tamanho);
                copiaMatriz(R, Rparcial, tamanho);
        } else{
                produtoMatrizes(Hvn, Rparcial, R, tamanho);
                copiaMatriz(R, Rparcial, tamanho);
        }
    }
}

// dada uma matriz A, a função calculaAutovalores_e_Autovetores produz a matriz de autovalores
// e de autovetores utilizando a fatoração QR para um dado epsilon

void calculaAutovalores_e_Autovetores (double A[MAX][MAX], double autovalores[MAX][MAX], double autovetores[MAX][MAX], 
                                        double epsilon, int tamanho){
    double produto[MAX][MAX], Q[MAX][MAX], R[MAX][MAX];
    double ehMaior = 1, temMaior;

    geraIdentidade(produto, tamanho);

    while (ehMaior == 1){

        fatoraQR (A, Q, R, tamanho); //fatora a matriz Ak
        produtoMatrizes (R, Q, A, tamanho); // gera a matriz Ak+1
        produtoMatrizes (produto, Q, autovetores, tamanho); //faz o produto das matrizes Q para chegar na matriz de autovetores
        copiaMatriz(autovetores, produto, tamanho);

        temMaior = 0;

        for (int i = 1; i < tamanho; i++){
            for (int j = 0; j < i; j++){

                if (fabs(A[i][j]) > epsilon) temMaior++;

            }
        }

        if (temMaior == 0) ehMaior = 0; //nao tem mais nenhum elemento abaixo da diagonal principal maior que epsilon
    }
    copiaMatriz(A, autovalores, tamanho);
}

//programa principal

int main() {

    double A[MAX][MAX], autovalores[MAX][MAX], autovetores[MAX][MAX], epsilon;
    int tamanho;

    lerMatriz(A, &tamanho);
    
    printf("Define the value of epsilon:\n");
    scanf("%lf", &epsilon);
    printf("\nCalculating eigenvalues ​​and eigenvectors...\n\n");

    calculaAutovalores_e_Autovetores(A, autovalores, autovetores, epsilon, tamanho);

    printf("Matrix of eigenvalues for epsilon = %lf :\n\n", epsilon);

    imprimeMatriz(autovalores, tamanho);

    printf("Matrix of eigenvectors for epsilon = %lf :\n\n", epsilon);

    imprimeMatriz(autovetores, tamanho);
    
    imprimeAutovaloresAutovetores(autovalores, autovetores, tamanho);

    printf("\nEnd of program.\n\n\n");
    
    return 0;
}
