import org.apache.commons.math3.special.Gamma;

import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

import static java.lang.Math.*;
import static org.apache.commons.math3.special.Erf.erf;
import static org.apache.commons.math3.special.Gamma.gamma;

public class Main {

    private static int N = 1000; // count of random variables
    private static int K = 128;  // MacLaren-Marsaglia Method box
    private static double KA = 1.35;  // 0.05
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


    private static double Hi2TestExp(double[] array, double a) {
        sort(array);
        int i = 0;
        int count, j;
        double xi2 = 0;
        double last = array[array.length-1];
        int intervalsCount = 10;
        double m = last/intervalsCount;

        for (j = 1; j <= intervalsCount; j++) {
            count = 0;
            double p = 0;
            while ((i < array.length) && array[i] < j*m ) {
                i++;
                count++;
            }
            p =  exp(-a*j*m);
            xi2 += Math.pow( (count - ((double) array.length) * p) , 2) /((double) array.length) * p;
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

    private static double Hi2TestForGeometric(int[] array, double p) {
        Map<Integer, Integer> count = new HashMap<>();


        for (int value : array) {
            if(count.containsKey(value)){
                count.put(value, count.get(value)+1);
            }else{
                count.put(value, 1);
            }
        }
        System.out.println("Freedom number: "+ count.size());
        double hi = 0;
        double d;
        for (int v : count.keySet()) {
            d = Math.pow(1-p, v-1)*p;
            hi += Math.pow(count.get(v) - N * d, 2) / N * d;
        }

        return hi;
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
        showIsAccepted("sqrt(n)Dn", resDn * sqrt(N), deltaKolm);

        double resDeltaHi = Hi2Test(Arrays.copyOf(alfas, alfas.length), 10);
        showIsAccepted("HI2", resDeltaHi, deltaHi);

        System.out.println(Arrays.toString(alfas));
        /*for (final double alfa : alfas) {
            System.out.println(alfa);
        }*/
        System.out.println();
    }

    private static void showResults(String what, double[] alfas, double realVariance, double realME) {
        System.out.println("**** " + what + "*****");

        System.out.println(Arrays.toString(alfas));
        int a = 0;
        for (final double alfa : alfas) {
            a += alfa;
        }
        double ueme = (double) a / N;
        System.out.println("unbiased estimate of mathematical expectation: " + ueme);
        System.out.println("real mathematical expectation: " + realME);

        double uev = 0.0;
        for (final double alfa : alfas) {
            uev += Math.pow(alfa - ueme, 2);
        }
        uev = uev / (N - 1);
        System.out.println("unbiased estimate of the variance: " + uev);
        System.out.println("real estimate of the variance: " + realVariance);

    }

    public static int[] generateBinomial(int size, double p, int m) {
        int[] result = new int[size];

        int aplha;
        int betta;

        for (int i = 0; i < size; ++i) {
            aplha = ThreadLocalRandom.current().nextInt(1000, 1000000);
            betta = ThreadLocalRandom.current().nextInt(1000, 1000000);
            aplha = aplha%2==0 ? aplha-1: aplha;
            betta = betta%2==0 ? betta-1: betta;
            result[i] = getBinomial(multiplicativeCongruentMethod(m, aplha, betta, M), p);
        }

        return result;
    }

    public static int[] generatePoison(int size, double lambda) {
        return generateBinomial(size, lambda / 10, 10);
    }

    public static int[] generateGeometric(int size, double p) {
        int[] result = new int[size];

        int aplha;
        int betta;

        for (int i = 0; i < size; ++i) {
            aplha = ThreadLocalRandom.current().nextInt(10000, 100000);
            betta = ThreadLocalRandom.current().nextInt(10000, 100000);
            aplha = aplha%2==0 ? aplha-1: aplha;
            betta = betta%2==0 ? betta-1: betta;
            result[i] = getGeometric(multiplicativeCongruentMethod(1, aplha, betta, M)[0], p);
        }

        return result;
    }

    public static int getGeometric(double value, double p) {
        return (int) Math.ceil(log(value)/ log(1-p));
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

    static void kolmogorovTestForNormal(double mu, double sg_2, double[] normdist)
    {
        double fotx, temp, d = 0.0, ka = sqrt(N);

        sort(normdist);

        for (int i = 0; i < N; i++)
        {
            fotx = 0.5 * (1 + erf((normdist[i] - mu) / (sqrt(2.0 * sg_2))));

            temp = abs((double)(i + 1) / N - fotx);

            if (d < temp)
                d = temp;
        }

        ka *= d;

        if (ka < KA)
            System.out.println( "Kolmogorov test: " + ka + " < " + KA + " => true sequence");
        else
            System.out.println( "Kolmogorov test: " + ka + " > " + KA + " => false sequence");
    }


    static double[] normalDistribution(double m, double s_2, int size)
    {
        double[] res = new double[size];
        int aplha;
        int betta;
        for (int i = 0; i < size; i++)
        {
            double temp = 0;
            for (int j = 0; j < 12; j++) {
                aplha = ThreadLocalRandom.current().nextInt(10000, 1000000);
                betta = ThreadLocalRandom.current().nextInt(10000, 1000000);
                aplha = aplha%2==0 ? aplha-1: aplha;
                betta = betta%2==0 ? betta-1: betta;
                temp += multiplicativeCongruentMethod(1, aplha, betta,M)[0];
            }
            temp -= 6;
            res[i] = m + sqrt(s_2) * temp;
        }

        return res;
    }


    static void kolmogorovTestForExponential(double a, double[] expdist)
    {
        double fotx, temp, d = 0.0, ka = sqrt(N);

        sort(expdist);

        for (int i = 0; i < N; i++)
        {
            fotx = 1.0 - exp(-a * expdist[i]);

            temp = abs((double)(i + 1) / N - fotx);

            if (d < temp)
                d = temp;
        }

        ka *= d;

        if (ka < KA)
            System.out.println("Kolmogorov test: " +ka +" < " +KA + " => true sequence");
        else
            System.out.println( "Kolmogorov test: " + ka + " > " + KA + " => false sequence" );
    }


    static double[] exponentialDistribution(double lm, int size)
    {
        double[] res = new double [size];
        int aplha;
        int betta;
        for (int i = 0; i < size; i++)
        {
            aplha = ThreadLocalRandom.current().nextInt(10000, 1000000);
            betta = ThreadLocalRandom.current().nextInt(10000, 1000000);
            aplha = aplha%2==0 ? aplha-1: aplha;
            betta = betta%2==0 ? betta-1: betta;
            res[i] = (-1.0 / lm) * log(multiplicativeCongruentMethod(1, aplha, betta,M)[0]);
        }

        return res;
    }

    static void kolmogorovTestForWeibull(double a, double b, double[] expdist)
    {
        double fotx, temp, d = 0.0, ka = sqrt(N);

        sort(expdist);

        for (int i = 0; i < N; i++)
        {
            fotx = 1.0 - exp(-a * pow(expdist[i], b));

            temp = abs((double)(i + 1) / N - fotx);

            if (d < temp)
                d = temp;
        }

        ka *= d;

        if (ka < KA)
            System.out.println("Kolmogorov test: " +ka +" < " +KA + " => true sequence");
        else
            System.out.println( "Kolmogorov test: " + ka + " > " + KA + " => false sequence" );
    }


    static double[] weibullDistribution(double a,double b,  int size)
    {
        double[] res = new double [size];
        int aplha;
        int betta;
        for (int i = 0; i < size; i++)
        {
            aplha = ThreadLocalRandom.current().nextInt(10000, 1000000);
            betta = ThreadLocalRandom.current().nextInt(10000, 1000000);
            aplha = aplha%2==0 ? aplha-1: aplha;
            betta = betta%2==0 ? betta-1: betta;
            res[i] = pow( (-1.0 / a) * log(multiplicativeCongruentMethod(1, aplha, betta,M)[0]), 1.0/b);
        }

        return res;
    }

    private static double Hi2TestWeibull(double[] array, double a) {
        sort(array);
        int i = 0;
        int count, j;
        double xi2 = 0;
        double last = array[array.length-1];
        int intervalsCount = 10;
        double m = last/intervalsCount;

        for (j = 1; j <= intervalsCount; j++) {
            count = 0;
            double p = 0;
            while ((i < array.length) && array[i] < j*m ) {
                i++;
                count++;
            }
            p =  exp(-a*j*m);
            xi2 += Math.pow( (count - ((double) array.length) * p) , 2) /((double) array.length) * p;
        }
        return xi2;
    }



    public static void main(String[] args) {
/*        double[] alfasKolm = multiplicativeCongruentMethod(N, paramAlfa, paramBetta, M);
        showResults("Multiplicative Congruent method", alfasKolm);

        double[] alfasMM = MacLarenMarsagliaMethod(N);
        showResults("MacLaren-Marsaglia method", alfasMM);*/

        final int m = 10;
        final double p = 0.05;
        final double lambda = 0.5;
        double[] valuesNormal = normalDistribution(4, 25, 1000);
        showResults("Normal method", valuesNormal, 25, 4);
        kolmogorovTestForNormal(4, 25, valuesNormal);

        double[] valuesExp = exponentialDistribution(0.5, 1000);
        showResults("Exp method", valuesExp, 1/0.5/0.5, 1/0.5);
        kolmogorovTestForExponential(0.5, valuesExp);
        System.out.println("Hi2 value: "+Hi2TestExp(valuesExp, 0.5));

        double[] valuesWei = weibullDistribution(4,0.5, 1000);
        showResults("Exp method", valuesWei, pow(4, -2/0.5)*(gamma(1 + 2 / 0.5) - pow(gamma(1+1/0.5), 2)), pow(4, -1/0.5)*gamma(1 + 1/ 0.5));
        kolmogorovTestForWeibull(4,0.5, valuesWei);

    }
}
