package com.guzova;

import java.io.File;
import java.util.Arrays;
import java.util.Scanner;

public class Main {

    public static double[][] mul(double[][] m1, double[][] m2) throws Exception {
        if (m1[0].length != m2.length)
            throw new Exception("Размеры матриц для перемножения не совпадают!");
        int max_k = (int) m1[0].length;
        double[][] result = new double[m1.length][m2[0].length];
        for (int i = 0; i < result.length; ++i) {
            for (int j = 0; j < result[0].length; ++j) {
                result[i][j] = 0;
                for (int k = 0; k < max_k; ++k) {
                    result[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        return result;
    }

    public static double[] mul(double[][] m, double[] v) throws Exception {
        if (m[0].length != v.length)
            throw new Exception("Размер матрицы и вектора для перемножения не совпадают!");
        int max_k = (int) v.length;
        double[] result = new double[m.length];
        for (int i = 0; i < result.length; ++i) {
            result[i] = 0;
            for (int k = 0; k < max_k; ++k) {
                result[i] += m[i][k] * v[k];
            }
        }
        return result;
    }

    public static double[][] sub(double[][] m1, double[][] m2) throws Exception {
        if (m1[0].length != m2[0].length || m1.length != m2.length)
            throw new Exception("Размеры матриц  не совпадают!");
        double[][] result = new double[m1.length][m1[0].length];
        for (int i = 0; i < result.length; ++i) {
            for (int j = 0; j < result[i].length; ++j) {
                result[i][j] = m1[i][j] - m2[i][j];
            }
        }
        return result;
    }

    public static double[] sub(double[] v1, double[] v2) throws Exception {
        if (v1.length != v2.length)
            throw new Exception("Размеры  векторов  не совпадают!");
        double[] result = new double[v1.length];
        for (int i = 0; i < result.length; ++i) {
            result[i] = v1[i] - v2[i];
        }
        return result;
    }
    //норма матрицы
    public static double norm(double[][] matrix) {
        double norm = 0;
        for (double[] row : matrix) {
            double summ = 0;
            for (double element : row)
                summ += Math.abs(element);
            if (summ > norm)
                norm = summ;
        }
        return norm;
    }

    public static double norm(double[] vector) {
        double norm = 0;
        for (double element : vector) {
            if (Math.abs(element) > norm)
                norm = Math.abs(element);
        }
        return norm;
    }

    public static void main(String[] args) throws Exception {
        File file = new File("D:\\proga_3sem_java\\vma1\\input.txt");
        Scanner input = new Scanner(file);
        int n = input.nextInt();
        double[][] A = new double[n][n];
        double[] b = new double[n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                A[i][j] = Double.parseDouble(input.next());
        for (int i = 0; i < n; ++i)
            b[i] = Double.parseDouble(input.next());
        double[][] E = new double[n][n];
        for (int d = 0; d < n; d++)
            E[d][d] = 1.0;

        System.out.println("Система: ");
        for (double[] row : A)
            System.out.println(Arrays.toString(row));
        System.out.println("Неоднородность: ");
        System.out.println(Arrays.toString(b));
        Solver output = Solver.solve(A, b);
        System.out.println("Определитель: " + output.det);
        System.out.println("Решение: ");
        System.out.println(Arrays.toString(output.x));
        System.out.println("Обратная матрица: ");
        for (double[] row : output.inv)
            System.out.println(Arrays.toString(row));
        System.out.println("Невязка: r = Ax-b");
        double[] r = sub(mul(A, output.x), b);
        System.out.println(Arrays.toString(r));
        System.out.println("Норма: ||r|| = " + norm(r));
        System.out.println("Невязка: R = A*A_inv-E");
        double[][] R = sub(mul(output.inv, A), E);
        for (double[] row : R)
            System.out.println(Arrays.toString(row));
        System.out.println("Норма: ||R|| = " + norm(R));
        System.out.println("Число обусловленности: " + norm(A) * norm(output.inv));
    }
}
