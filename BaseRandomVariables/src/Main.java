import java.util.Arrays;

public class Main {

    private static int N = 1000;
    private static int K = 128;
    private static int paramBetta = 131075;
    private static int paramAlfa = 131075;
    private static double M = Math.pow(2, 10);

    public static double mod(double what, double module) {
        return what - module * (Math.floor( what / module));
    }

    public static double[] MacLarenMarsagliaMethod(int size) {
        double result[] = new double[size];

        int BRVsize = size*K;
        double[] b = multiplicativeCongruentMethod(BRVsize, paramAlfa, paramBetta, M);
        double[] c = multiplicativeCongruentMethod(BRVsize, 79507, 79507, M);
        double[] V = new double[K];

        for (int i = 0; i < size; i++) {
            System.arraycopy(b, (K*i), V, 0, K);

            result[i] = V[(int) Math.floor(c[i] * K)];
        }

        return result;
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

    public static void main(String[] args) {
        double[] alfas = MacLarenMarsagliaMethod(N);

        System.out.println(Arrays.toString(alfas));
    }
}
