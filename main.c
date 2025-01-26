#include <stdio.h>
#include <math.h>

#define TABLE_SIZE 28

// Tabelas fornecidas no código MATLAB
const double angles[TABLE_SIZE] = {
    0.78539816339745, 0.46364760900081, 0.24497866312686, 0.12435499454676,
    0.06241880999596, 0.03123983343027, 0.01562372862048, 0.00781234106010,
    0.00390623013197, 0.00195312251648, 0.00097656218956, 0.00048828121119,
    0.00024414062015, 0.00012207031189, 0.00006103515617, 0.00003051757812,
    0.00001525878906, 0.00000762939453, 0.00000381469727, 0.00000190734863,
    0.00000095367432, 0.00000047683716, 0.00000023841858, 0.00000011920929,
    0.00000005960464, 0.00000002980232, 0.00000001490116, 0.00000000745058
};

const double Kvalues[TABLE_SIZE] = {
    0.70710678118655, 0.63245553203368, 0.61357199107790, 0.60883391251775,
    0.60764825625617, 0.60735177014130, 0.60727764409353, 0.60725911229889,
    0.60725447933256, 0.60725332108988, 0.60725303152913, 0.60725295913894,
    0.60725294104140, 0.60725293651701, 0.60725293538591, 0.60725293510314,
    0.60725293503245, 0.60725293501477, 0.60725293501035, 0.60725293500925,
    0.60725293500897, 0.60725293500890, 0.60725293500889, 0.60725293500888
};

// Função CORDIC
void cordic(double beta, int n, double *cosine, double *sine, double *tangent) {
    // Ajustar beta para o intervalo principal [-pi/2, pi/2]
    int flip_sign = 0;
    if (beta < -M_PI_2 || beta > M_PI_2) {
        if (beta < 0)
            beta += M_PI;
        else
            beta -= M_PI;
        flip_sign = 1;  // O sinal precisa ser invertido
    }

    // Inicialização de variáveis
    double x = 1.0, y = 0.0;  // Vetor inicial [1, 0]
    double power_of_two = 1.0;
    double angle = angles[0];

    // Iteração do método CORDIC
    for (int i = 0; i < n; i++) {
        int sigma = (beta < 0) ? -1 : 1;
        double factor = sigma * power_of_two;

        // Rotação
        double x_new = x - factor * y;
        double y_new = y + factor * x;
        x = x_new;
        y = y_new;

        // Atualização do ângulo
        beta -= sigma * angle;
        power_of_two /= 2;

        // Atualizar o próximo ângulo
        if (i + 1 < TABLE_SIZE)
            angle = angles[i + 1];
        else
            angle /= 2;
    }

    // Aplicar o fator de escala
    double Kn = (n < TABLE_SIZE) ? Kvalues[n - 1] : Kvalues[TABLE_SIZE - 1];
    x *= Kn;
    y *= Kn;

    // Ajustar os valores de saída
    if (flip_sign) {
        x = -x;
        y = -y;
    }

    *cosine = x;
    *sine = y;

    // Cálculo da tangente
    *tangent = y / x;
}

// Função principal para testar
int main() {
    double beta = 3.141592619 / 6;  // Ângulo em radianos
    int n = 20;         // Número de iterações
    double cosine, sine, tangent;

    cordic(beta, n, &cosine, &sine, &tangent);

    printf("Cosseno: %.12f\n", cosine);
    printf("Seno: %.12f\n", sine);
    printf("Tangente: %.12f\n", tangent);

    return 0;
}
