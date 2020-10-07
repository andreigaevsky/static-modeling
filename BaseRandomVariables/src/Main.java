import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;

public class Main {

    private static int N = 1000; // count of random variables
    private static int K = 128;  // MacLaren-Marsaglia Method box
    private static double deltaKolm = 1.36; // Kolmogorov distribution function value
    private static double deltaHi = 16.92; // chi square distribution function value
    private static double deltaHiForBernoulli = 3.84; // chi square distribution function value
    private static int paramBetta = 131075; // any value
    private static int paramAlfa = 131075; // any value
    private static double M = Math.pow(2, 31);

    private static void sort(double[] array) {
        double temp;
        int i, j;
        for (i = 0; i < array.length; i++) {
            for (j = 1; j < array.length - i; j++) {
                if (array[j - 1] > array[j]) {
                    temp = array[j];
                    array[j] = array[j - 1];
                    array[j - 1] = temp;
                }
            }
        }
    }

    private static double KolmogorovTest(double[] array) {
        sort(array);
        int i = 0;
        double D = 0;
        double temp;
        for (i = 0; i < array.length; i++) {
            temp = Math.abs(((double) (i + 1)) / array.length - array[i]);
            if (D < temp) {
                D = temp;
            }
        }
        return D;
    }

    private static double Hi2Test(double[] array, int intervalsCount) {
        sort(array);
        int i = 0;
        int count, j;
        double xi2 = 0;
        for (j = 1; j <= intervalsCount; j++) {
            count = 0;
            while ((i < array.length) && array[i] < (double) j / intervalsCount) {
                i++;
                count++;
            }
            xi2 += Math.pow((count - ((double) array.length) / intervalsCount), 2) / (((double) array.length) / intervalsCount);
        }
        return xi2;
    }

    private static double Hi2TestForBernoulli(int[] array, double p) {
        int zero = 0;
        int one = 0;
        for (int value : array) {
            if (value == 0) {
                ++zero;
            } else {
                one++;
            }
        }
        double hiFirst = Math.pow(zero - N * (1 - p), 2) / N * (1 - p);
        double hiSecond = Math.pow(one - N * p, 2) / N * p;
        return hiFirst + hiSecond;
    }

    private static double Hi2TestForBinomial(int[] array, double p, int maxValue) {
        int[] count = new int[maxValue+1];
        Arrays.fill(count, 0);
        for (int value : array) {
            count[value] = ++count[value];
        }
        double hi = 0;
        double v;
        for (int i = 0; i <= maxValue; i++) {
            v = getBinomialValue(i, maxValue, p);
            hi += Math.pow(count[i] - N * v, 2) / N * v;
        }

        return hi;
    }

    private static double getBinomialValue(int value, int m, double p){
        return (double)factorial(m)/factorial(m-value)/factorial(value)*Math.pow(p, value)*Math.pow(1-p,m-value);
    }

    public static long factorial(int number) {
        long result = 1;

        for (int factor = 2; factor <= number; factor++) {
            result *= factor;
        }

        return result;
    }


    public static double mod(double what, double module) {
        return what - module * (Math.floor(what / module));
    }

    public static double[] MacLarenMarsagliaMethod(int size) {
        double result[] = new double[size];

        int BRVsize = size + K;
        double[] b = multiplicativeCongruentMethod(BRVsize, paramAlfa, paramBetta, M);
        double[] c = multiplicativeCongruentMethod(size, 79507, 79507, M);
        double[] V = new double[K];
        System.arraycopy(b, 0, V, 0, K);
        int s;
        for (int i = 0; i < size; i++) {

            s = (int) Math.floor(c[i] * K);
            result[i] = V[s];
            V[s] = b[i + K];
        }

        return result;
    }

    private static void showIsAccepted(String what, double result, double expected) {
        System.out.println(what + " = " + result + " < " + expected + " is " + (result < expected));
    }

    public static double[] multiplicativeCongruentMethod(int size, double paramAlfa, double paramBetta, double M) {
        double[] result = new double[size];

        double nextParamAlfa = mod(paramAlfa * paramBetta, M);
        for (int i = 0; i < size; i++) {
            result[i] = nextParamAlfa / M;
            nextParamAlfa = mod(nextParamAlfa * paramBetta, M);
        }
        return result;
    }

    private static void showResults(String what, double[] alfas) {
        System.out.println("**** " + what + "*****");

        double resDn = KolmogorovTest(Arrays.copyOf(alfas, alfas.length));
        showIsAccepted("sqrt(n)Dn", resDn * Math.sqrt(N), deltaKolm);

        double resDeltaHi = Hi2Test(Arrays.copyOf(alfas, alfas.length), 10);
        showIsAccepted("HI2", resDeltaHi, deltaHi);

        System.out.println(Arrays.toString(alfas));
        /*for (final double alfa : alfas) {
            System.out.println(alfa);
        }*/
        System.out.println();
    }

    private static void showResults(String what, int[] alfas, double realVariance, double realME) {
        System.out.println("**** " + what + "*****");

        System.out.println(Arrays.toString(alfas));
        int a = 0;
        for (final int alfa : alfas) {
            a += alfa;
        }
        double ueme = (double) a / N;
        System.out.println("unbiased estimate of mathematical expectation: " + ueme);
        System.out.println("real mathematical expectation: " + realME);

        double uev = 0.0;
        for (final int alfa : alfas) {
            uev += Math.pow(alfa - ueme, 2);
        }
        uev = uev / (N - 1);
        System.out.println("unbiased estimate of the variance: " + uev);
        System.out.println("real estimate of the variance: " + realVariance);

        System.out.println();
    }

    public static int[] generateBinomial(int size, double p, int m) {
        int[] result = new int[size];

        int aplha;
        int betta;

        for (int i = 0; i < size; ++i) {
            aplha = ThreadLocalRandom.current().nextInt(10000, 100000);
            betta = ThreadLocalRandom.current().nextInt(10000, 100000);
            aplha = aplha%2==0 ? aplha-1: aplha;
            betta = betta%2==0 ? betta-1: betta;
            result[i] = getBinomial(multiplicativeCongruentMethod(m, aplha, betta, M), p);
        }

        return result;
    }

    public static int[] generatePoison(int size, double lambda) {
        return generateBinomial(size, lambda / 10, 10);
    }

    public static int[] generateBernoulli(int size, double p) {
        int[] result = new int[size];

        int aplha;
        int betta;

        for (int i = 0; i < size; ++i) {
            aplha = ThreadLocalRandom.current().nextInt(10000, 100000);
            betta = ThreadLocalRandom.current().nextInt(10000, 100000);
            aplha = aplha%2==0 ? aplha-1: aplha;
            betta = betta%2==0 ? betta-1: betta;
            result[i] = getBernoulli(multiplicativeCongruentMethod(1, aplha, betta, M)[0], p);
        }

        return result;
    }

    public static int getBernoulli(double value, double p) {
        return value > p ? 0 : 1;
    }

    public static int getBinomial(double[] array, double p) {
        int x = 0;
        for (double v : array) {
            if (p > v) {
                ++x;
            }
        }
        return x;
    }

    public static void main(String[] args) {
/*        double[] alfasKolm = multiplicativeCongruentMethod(N, paramAlfa, paramBetta, M);
        showResults("Multiplicative Congruent method", alfasKolm);

        double[] alfasMM = MacLarenMarsagliaMethod(N);
        showResults("MacLaren-Marsaglia method", alfasMM);*/

        final int m = 10;
        final double p = 0.05;
        final double lambda = 0.5;
        int[] valuesPoison = generateBinomial(N,p,m);
        showResults("Poison method", valuesPoison, m*p * (1 - p), m*p);
        double resDeltaHiPoison = Hi2TestForBinomial(Arrays.copyOf(valuesPoison, valuesPoison.length), p, m);
        showIsAccepted("HI2", resDeltaHiPoison, deltaHi);

        int[] valuesBernoulli = generateBernoulli(N, 0.6);
        showResults("Bernoulli method", valuesBernoulli, 0.6 * (1 - 0.6), 0.6);
        double resDeltaHi = Hi2TestForBernoulli(Arrays.copyOf(valuesBernoulli, valuesBernoulli.length), 0.6);
        showIsAccepted("HI2", resDeltaHi, deltaHiForBernoulli);
    }
}
