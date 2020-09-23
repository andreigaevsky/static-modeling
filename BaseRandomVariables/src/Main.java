import java.util.Arrays;

public class Main {

    private static int N = 1000; // count of random variables
    private static int K = 128;  // MacLaren-Marsaglia Method box
    private static double deltaKolm = 1.36; // Kolmogorov distribution function value
    private static double deltaHi = 16.92; // chi square distribution function value
    private static int paramBetta = 131075; // any value
    private static int paramAlfa = 131075; // any value
    private static double M = Math.pow(2, 31);

    private static void sort(double[] array) {
        double temp;
        int i, j;
        for (i = 0; i < array.length; i++)
            for (j = 1; j < array.length - i; j++) {
                if (array[j - 1] > array[j]) {
                    temp = array[j];
                    array[j] = array[j - 1];
                    array[j - 1] = temp;
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
            if (D < temp)
                D = temp;
        }
        return D;
    }


    private static double Hi2Test(double[] array, int intervalsCount) {
         sort(array);
        int i = 0;
        int count = 0, j = 0;
        double xi2 = 0;
        for (j = 1; j < intervalsCount; j++) {
            count = 0;
            while ((i < array.length) && (array[i] < ((double) j) / intervalsCount)) {
                i++;
                count++;
            }
            xi2 += Math.pow((count - ((double)array.length) / intervalsCount), 2) / (((double)array.length) / intervalsCount);
        }
        return xi2;
    }

    public static double mod(double what, double module) {
        return what - module * (Math.floor( what / module));
    }

    public static double[] MacLarenMarsagliaMethod(int size) {
        double result[] = new double[size];

        int BRVsize = size+K;
        double[] b = multiplicativeCongruentMethod(BRVsize, paramAlfa, paramBetta, M);
        double[] c = multiplicativeCongruentMethod(size, 79507, 79507, M);
        double[] V = new double[K];
        System.arraycopy(b, 0, V, 0, K);
        int s;
        for (int i = 0; i < size; i++) {

            s = (int) Math.floor(c[i] * K);
            result[i] = V[s];
            V[s] = b[i+K];
        }

        return result;
    }

    private static void showIsAccepted(String what, double result, double expected) {
        System.out.println(what+" = " + result + " < " +expected+" is "+ (result < expected));
    }

    public static double[] multiplicativeCongruentMethod(int size, int paramAlfa, int paramBetta, double M) {
        double[] result = new double[size];

        double nextParamAlfa = mod(paramAlfa * paramBetta, M);
        for (int i = 0; i < size; i++) {
            result[i] = nextParamAlfa / M;
            nextParamAlfa = mod(nextParamAlfa * paramBetta, M);
        }
        return result;
    }

    private static void showResults (String what, double [] alfas) {
        System.out.println("**** "+what+"*****");

        double resDn = KolmogorovTest(Arrays.copyOf(alfas, alfas.length));
        showIsAccepted("sqrt(n)Dn", resDn*Math.sqrt(N)*resDn, deltaKolm);

        double resDeltaHi = Hi2Test(Arrays.copyOf(alfas, alfas.length),10);
        showIsAccepted("HI2", resDeltaHi, deltaHi);

        System.out.println(Arrays.toString(alfas));
        /*for (final double alfa : alfas) {
            System.out.println(alfa);
        }*/
        System.out.println();
    }

    public static void main(String[] args) {
        double[] alfasKolm = multiplicativeCongruentMethod(N, paramAlfa, paramBetta, M);
        showResults("Multiplicative Congruent method", alfasKolm);

        double[] alfasMM = MacLarenMarsagliaMethod(N);
        showResults("MacLaren-Marsaglia method", alfasMM);
    }
}