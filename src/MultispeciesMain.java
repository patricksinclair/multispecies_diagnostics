public class MultispeciesMain {
    public static void main(String[] args) {

        //todo - make sure alpha is changed to correspond with the new parameter values when relevant
        int nReps = Integer.parseInt(args[0]);
        double tau_val = 0.2;

        //these ints tell us which chunks of parameter values to run
        //possible values range from 1 to 20.
        int K_startIndex = 1, K_endIndex = 4; //chunk 1 values
        //int K_startIndex = 5, K_endIndex = 8; //chunk 2 values
        //int K_startIndex = 9, K_endIndex = 12; //chunk 3 values
        //int K_startIndex = 13, K_endIndex = 16; //chunk 4 values
        //int K_startIndex = 17, K_endIndex = 20; //chunk 5 values


        BioSystem.varyingDeteriorationAndThreshold(nReps, tau_val, K_startIndex, K_endIndex);
    }
}
