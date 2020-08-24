public class MultispeciesMain {
    public static void main(String[] args) {

        //todo - make sure alpha is changed to correspond with the new parameter values when relevant
        int nReps = Integer.parseInt(args[0]);
        double tau_val = 0.2;

        //these ints tell us which chunks of parameter values to run
        //possible values range from 0 to 19.
        //int K_startIndex = 0, K_endIndex = 3; //chunk 1 values
        //something went wrong with chunks 2 and 3, so changed the indexes there to reflect that
        //int K_startIndex = 6, K_endIndex = 7; //chunk 2 values
        //int K_startIndex = 10, K_endIndex = 11; //chunk 3 values
        //int K_startIndex = 12, K_endIndex = 15; //chunk 4 values
        //int K_startIndex = 16, K_endIndex = 19; //chunk 5 values

        //BioSystem.varyingDeteriorationAndThreshold(nReps, tau_val, K_startIndex, K_endIndex);
        BioSystem.varyingDeteriorationAndThresholdN_extraDetRatios(nReps, tau_val);
    }
}
