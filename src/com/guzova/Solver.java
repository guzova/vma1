package com.guzova;

import java.util.Arrays;

public class Solver {
    private Solver(int n) {
        this.inv = new double[n][n];
        this.x = new double[n];
        this.det = 1;
    }

    public double det;
    public double[][] inv;
    public double[] x;

    public static Solver solve(double[][] _A, double[] _b) throws Exception {
        int n = _A.length;
        double[][] A = new double[n][n];
        double[] b = new double[n];
        System.arraycopy(_b, 0, b, 0, n);
        for (int d = 0; d < n; d++)
            System.arraycopy(_A[d], 0, A[d], 0, n);
        for (double[] row : A)
            if (row.length != n)
                throw new Exception("Неправильный размер");
        if (b.length != n)
            throw new Exception("Неправильный размер ");

        Solver output = new Solver(n);
        double[][] E = new double[n][n];
        for (int d = 0; d < n; d++)
            E[d][d] = 1.0;

        for (int k = 0; k < n; k++) {

            //в столбце находится максимальный по модулю элемент (ведущий)
            int max_ind = k;
            double max_value = 0;
            for (int i = k; i < n; i++) {
                if (Math.abs(A[i][k]) > Math.abs(max_value)) {
                    max_ind = i;
                    max_value = A[i][k];
                }
            }
            // и вместе со строкой переставляется наверх
            if (max_ind != k) {

                double[] temp = A[k];
                A[k] = A[max_ind];
                A[max_ind] = temp;

                double[] _temp = E[k];
                E[k] = E[max_ind];
                E[max_ind] = _temp;

                double __temp = b[k];
                b[k] = b[max_ind];
                b[max_ind] = __temp;

                output.det *= -1;
            }
            if (Math.abs(A[k][k]) < 1e-15)
                throw new Exception("Определитель около 0!");
            output.det *= A[k][k];

            //делим верхнюю строку, неоднородность, строку обратной матрицы на ведущий эдемент
            for (int i = k + 1; i < n; i++) {
                A[k][i] /= A[k][k];
            }
            b[k] /= A[k][k];
            for (int i = 0; i < n; i++) {
                E[k][i] /= A[k][k];
            }
            A[k][k] = 1;

            //отнимаем ведущую строку от нижних домноженную на первые элементы этих строк
            for (int i = k + 1; i < n; i++) {
                for (int j = k + 1; j < n; j++) {
                    A[i][j] -= A[i][k] * A[k][j];
                }
                b[i] -= b[k] * A[i][k];
                for (int j = 0; j < n; j++) {
                    E[i][j] -= A[i][k] * E[k][j];
                }
                A[i][k] = 0;
            }
        }

        //обратный ход
        for (int k = n - 1; k >= 0; --k) {
            output.x[k] = b[k];
            System.arraycopy(E, 0, output.inv, 0, n);
            for (int i = n - 1; i > k; --i) {
                output.x[k] -= A[k][i] * output.x[i];
                for (int j = 0; j < n; ++j) {
                    output.inv[k][j] -= A[k][i] * output.inv[i][j];
                }
            }
        }
        return output;
    }
}

