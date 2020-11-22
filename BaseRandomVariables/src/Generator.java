import static java.lang.Math.log;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.util.concurrent.ThreadLocalRandom;

public class Generator {

    public static double[] multiplicativeCongruentMethod(int size, double paramAlfa, double paramBetta, double M) {
        double[] result = new double[size];

        double nextParamAlfa = mod(paramAlfa * paramBetta, M);
        for (int i = 0; i < size; i++) {
            result[i] = nextParamAlfa / M;
            nextParamAlfa = mod(nextParamAlfa * paramBetta, M);
        }
        return result;
    }

    private static double mod(double what, double module) {
        return what - module * (Math.floor(what / module));
    }

    public static int[] generateBinomial(int size, double p, int m, double M) {
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

    public static int[] generatePoison(int size, double lambda, double M) {
        return generateBinomial(size, lambda / 10, 10, M);
    }

    public static int[] generateGeometric(int size, double p, double M) {
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

    private static int getGeometric(double value, double p) {
        return (int) Math.ceil(log(value)/ log(1-p));
    }

    public static int[] generateBernoulli(int size, double p, double M) {
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

    private static int getBernoulli(double value, double p) {
        return value > p ? 0 : 1;
    }

    private static int getBinomial(double[] array, double p) {
        int x = 0;
        for (double v : array) {
            if (p > v) {
                ++x;
            }
        }
        return x;
    }

    public static double[] normalDistribution(double m, double s_2, int size, double M)
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

    public static double[] evenDistribution(double a, double b, int size, double M)
    {
        double[] res = new double[size];
        int aplha;
        int betta;
        double n;
        for (int i = 0; i < size; i++)
        {
            aplha = ThreadLocalRandom.current().nextInt(10000, 1000000);
            betta = ThreadLocalRandom.current().nextInt(10000, 1000000);
            aplha = aplha%2==0 ? aplha-1: aplha;
            betta = betta%2==0 ? betta-1: betta;
            n = multiplicativeCongruentMethod(1, aplha, betta,M)[0];
            res[i] = (b-a)*n + a;
        }

        return res;
    }

    public static double[] exponentialDistribution(double lm, int size, double M)
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

    public static double[] weibullDistribution(double a,double b,  int size, double M)
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
}
