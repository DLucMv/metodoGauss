/* Método de Gauss para a soulução de Sistemas Lineares.
 * Autor: Davison Lucas Mendes Viana - Engenharia da Computação - IFCE
 * Atividade desenvolvida para a disciplina de Cálculo Numérico
 * Professor: Glauber Ferreira Cintra */

#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <conio.h>
#include <windows.h>

enum DOS_COLORS {
    BLACK,
    BLUE,
    GREEN,
    CYAN,
    RED,
    MAGENTA,
    BROWN,
    LIGHT_GRAY,
    DARK_GRAY,
    LIGHT_BLUE,
    LIGHT_GREEN,
    LIGHT_CYAN,
    LIGHT_RED,
    LIGHT_MAGENTA,
    YELLOW,
    WHITE };

double **criaMatriz(int linha, int coluna);
void leMatriz(double **matriz, int linha, int coluna);
void exibeMatriz(double **matriz, int linha, int coluna);
void algoritmoGauss(double **m, int linha);
int sRetro(double **m, int linha, double x[]);
void type(int tipo, int linha, double x[]);
double *criarVetor(int linha);
void textcolor(int color);


int main() {

    setlocale(LC_ALL, "");

    char op = 'n';
    int n, tipo;
    double *x;
    double **mat;

    do {
        printf("\nQual a Dimens\xE3o da matriz aumentada (n x n + 1): \n");
        scanf("%d", &n);

        x = criarVetor(n);

        mat = criaMatriz(n, n + 1);
        if(mat == NULL){
            printf("Deu pau!\nMem\xF3ria insuficiente!\n");
            return 1;
        }

        leMatriz(mat, n, n + 1);

        printf("\n\tExibindo a Matriz (%dx%d)\n", n, n + 1);
        exibeMatriz(mat, n, n + 1);

        algoritmoGauss(mat, n);

        printf("\n\tExibindo a Matriz Normalizada (%dx%d)\n", n, n + 1);
        exibeMatriz(mat, n, n + 1);

        tipo = sRetro(mat, n, x);
        type(tipo, n, x);

        printf("\nDeseja verificar uma nova matriz? (s/n)\n");
        scanf("%s", &op);

    }while (op == 's');

    getch();
    return 0;
}


double **criaMatriz(int l, int c) {
    /* Se houver memória suficiente, cria uma matriz de double com l linhas *
     * e c colunas e devolve um ponteiro para a matriz; Caso contrário,     *
     * devolve um ponteiro nulo.                                            */

    int i, j;
    double **m = NULL;

    m = malloc(sizeof (double *) * l);

    textcolor(LIGHT_RED);
    if(m == NULL){
        printf("\nMem\xF3rio insuficiente.\n");
        return NULL;
    }

    for (i = 0; i < l; ++i) {
        m[i] = malloc(sizeof (double) * c);
        if(m[i] == NULL){
            printf("\nMem\xF3ria insuficiente.\n");
            for (j = 0; j < i; ++j) {
                free(m[j]);
            }
            free(m);
            return NULL;
        }
    }
    textcolor(WHITE);

    return m;

}


void leMatriz(double **m, int l, int c) {
    /* Lê valores para a matriz M, que é uma matriz double com l linhas e   *
     * c colunas alocada dinamicamente.                                     */

    int i, j;

    printf("\nInsira os elementos da matriz aumentada...\n");

    for (i = 0; i < l; ++i) {
        for (j = 0; j < c; ++j) {
            printf("M[%d][%d]: ", i + 1, j + 1);
            scanf("%lf", &m[i][j]);
        }
    }

    printf("\n");

}


void exibeMatriz(double **m, int l, int c) {
    /* Exibe o conteúdo de uma matriz de double com l linhas e c colunas    *
     * alocada dinamicamente.                                               */

    textcolor(LIGHT_BLUE);
    for (int i = 0; i < l; ++i) {
        printf("|");
        for (int j = 0; j < c; ++j) {
            printf("%10.3lf", m[i][j]);
        }
        printf("\t|\n");
    }
    textcolor(WHITE);

    printf("\n");
}


void algoritmoGauss(double **m, int l){
    /* Recebe m, a matriz aumentada de um SL com n variáveis e      *
     * transforma m na matriz aumentada de um SL TS quivalente ao   *
     * SL fornecido na entrada.                                     */

    int i, j, k;
    double mj;
    double *aux;
    double **mat = m;

    for (i = 0; i < l - 1; ++i) {
        if(mat[i][i] == 0){ /* Pivô Nulo. */
            j = i + 1;
            while (j < l && mat[j][i] == 0){
                j++;
            }
            if(j < l){
                aux = mat[i];
                mat[i] = mat[j];
                mat[j] = aux;
            }
        }
        if(mat[i][i] != 0){
            for(j = i + 1; j < l ; j++){
                mj = (-mat[j][i]) / (mat[i][i]);
                mat[j][i] = 0;
                for(k = i + 1; k < l + 1; k++){
                    mat[j][k] += (mj * mat[i][k]);
                }
            }
        }
    }

}


int sRetro(double **m, int l, double x[]){
    /* Recebe M, a matriz aumentada de um SL TS com l variáveis.            *
     * Se o SL for determindado, coloca em x a solução do SL e devolve 0;   *
     * Se for indeterminado, coloca em x uma solução e devolve 1;           *
     * Se for incompativél, devolve 2.                                      */

    int tipo = 0;
    double soma;

    for (int i = l - 1; i >= 0; i--) {
        soma = 0;
        for(int j = i + 1; j < l; j++){
            soma += m[i][j] * x[j];
        }
        if(m[i][i] == 0){
            if(m[i][l] == soma){ /* x[i] é variável livre. */
                x[i] = 0;
                tipo = 1;
            } else {
                return 2; /* SL incompatível. */
            }
        } else {
            x[i] = (m[i][l] - soma) / (m[i][i]);
        }
    }
    return tipo;

}


void type(int tipo, int l, double x[]){

    if(tipo == 2){
        printf("SL INCOMPAT\xCDVEL\n");
    } else {
        printf("SL %sDETERMINADO\n", tipo? "IN" : "");
        for(int i = 0; i < l; i++){
            printf("x[%d] = %10.3lf\n", i + 1, x[i]);
        }
    }

}


double *criarVetor(int l) {
    /* Se houver memória suficiente, cria um vetor de doubles com   *
     * tamanho l. Caso contrário, devolve um ponteiro nulo.         */

    double *x;
    x = malloc(sizeof(double) * l);

    textcolor(LIGHT_RED);
    if (x == NULL) {
        printf("Não foi possivel criar o vetor!\n");
        return NULL;
    }
    textcolor(WHITE);

    return x;

}


void textcolor (int color)
{
    static int __BACKGROUND;

    HANDLE h = GetStdHandle ( STD_OUTPUT_HANDLE );
    CONSOLE_SCREEN_BUFFER_INFO csbiInfo;


    GetConsoleScreenBufferInfo(h, &csbiInfo);

    SetConsoleTextAttribute (GetStdHandle (STD_OUTPUT_HANDLE),
                             color + (__BACKGROUND << 4));
}
