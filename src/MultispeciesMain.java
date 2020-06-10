public class MultispeciesMain {
    public static void main(String[] args) {

        int nReps = Integer.parseInt(args[0]);
        double tau_val = 0.2;
        BioSystem.varyingDeteriorationAndThreshold(nReps, tau_val);
    }
}
