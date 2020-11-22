import java.util.Arrays;
import java.util.Objects;

public class SystemResolver {
    private double[][] A;
    private double[] f;
    private double[] solution;
    private static final double M = Math.pow(2, 31);
    private int N;
    private int m;

    SystemResolver(double[][] A, double[] f, int N, int m) {
        validateSizes(A);
        this.A = shiftToLeft(A);
        this.f = Arrays.copyOf(f, f.length);
        this.N = N;
        this.m = m;
    }

    private double[][] shiftToLeft(double[][] matrix) {
        double[][] result = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; ++i){
            for(int j = 0; j < matrix[0].length; ++j){
                if(i ==  j){
                    result[i][j] = 1 + -1*matrix[i][j];
                }
                else{
                    result[i][j] = -1*matrix[i][j];
                }
            }
        }
        return result;
    }

    public void resolve() {
        this.solution = new double[A.length];

        double[][] P = fillP(A);
        double[][] H = getIdentityMatrix(f.length);

        double[] pi = new double[A.length];
        Arrays.fill(pi, 1d / A.length);

        int[] I = new int[N + 1];
        double[] Q = new double[N + 1];
        double[] ksi = new double[m];

        for (int solution = 0; solution < this.solution.length; ++solution) {
            Arrays.fill(ksi, 0d);
            Arrays.fill(I, 0);
            Arrays.fill(Q, 0);
            double[] h = H[solution];
            for (int j = 0; j < m; ++j) {
                for (int k = 0; k <= N; ++k) {
                    double alpha = Generator.evenDistribution(0, 1, 1, M)[0];
                    double cumSum = 0;
                    for (int c = 0; c < pi.length; ++c) {
                        cumSum += pi[c];
                        if (alpha < cumSum) {
                            I[k] = c;
                            break;
                        }
                    }
                }
                Q[0] = pi[I[0]] > 0 ? h[I[0]] / pi[I[0]] : 0;
                for (int k = 1; k <= N; ++k) {
                    Q[k] = 0 < P[I[k - 1]][I[k]] ? Q[k - 1] * A[I[k - 1]][I[k]] / P[I[k - 1]][I[k]] : 0;
                }
                for (int k = 0; k <= N; ++k)
                    ksi[j] = ksi[j] + Q[k] * f[I[k]];
            }
            this.solution[solution] = Arrays.stream(ksi).sum() /  ksi.length;
        }
    }

    public void showSolution() {
        displayArray(solution);
    }

    public double[] getDiscrepancy(double[] actualSolution) {
        double[] discrepancy = new double[solution.length];
        for (int i = 0; i < solution.length; ++i) {
            discrepancy[i] = actualSolution[i] - solution[i];
        }
        return discrepancy;
    }

    private void displayArray(double[] array) {
        Arrays.stream(array).forEach(System.out::println);
    }

    private double[][] fillP(double[][] A) {
        double[][] P = new double[A.length][A.length];
        for (int i = 0; i < A.length; ++i) {
            for (int j = 0; j < A.length; ++j) {
                if (A[i][j] == 0) {
                    P[i][j] = 0;
                } else {
                    P[i][j] = 1d / A.length;
                }
            }
        }
        return P;
    }

    private double [][] getIdentityMatrix (int size) {
        double [][] matrix = new double[size][size];
        for(int i = 0; i < size; ++i) {
            for(int j = 0; j < size; j++) {
                if(i == j) {
                    matrix[i][j] = 1;
                }else{
                    matrix[i][j] = 0;
                }
            }
        }
        return matrix;
    }

    private void validateSizes(double[][] matrix) {
        if (Objects.isNull(matrix) || matrix[0].length != matrix.length) {
            throw new IllegalArgumentException();
        }
    }
}
