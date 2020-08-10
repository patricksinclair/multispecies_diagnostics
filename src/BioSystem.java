import org.apache.commons.math3.distribution.PoissonDistribution;
import java.util.ArrayList;
import java.util.Random;
import java.util.stream.IntStream;

class BioSystem {

    private Random rand = new Random();
    //need to initialise these in the performAction method in order to incorporate the tau halving
    private PoissonDistribution poiss_imigration;
    private PoissonDistribution poiss_deterioration;
    private PoissonDistribution poiss_migration;
    private PoissonDistribution poiss_migration_edge;


    private double alpha, c_max; //steepness and max val of antimicrobial concn
    private double scale, sigma; //mic distb shape parameters
    private ArrayList<Microhabitat> microhabitats;
    private double time_elapsed, exit_time; //exit time is the time it took for the biofilm to reach the thickness limit, if it did
    private int immigration_index;

    private double deterioration_rate;
    private double biofilm_threshold;

    private int K = 550; //carrying capacity
    private double max_gRate = 0.083; //max growth rate =  2/day
    private double immigration_rate = 20.;
    private double migration_rate = 1.;
    private double tau;
    private double delta_z = 1.;
    private int thickness_limit = 22; //this is how big the system can get before we exit. should reduce overall simulation duration todo-change back to "50" for big runs
    private int detachments_counter = 0, deaths_counter = 0, replications_counter = 0, immigrations_counter = 0, tau_halves_counter = 0; //last one is the number of times tau had to be halved due to double events


    public BioSystem(double deterioration_ratio, double biofilm_threshold, double tau){

        //this constructor is used purely for the detachment rate determination in the biocide free environment
        this.alpha = 0.;
        this.c_max = 0.;
        //this scale and sigma correspond to 99% susceptible
        this.scale = 2.71760274;
        this.sigma = 0.56002833;
        this.microhabitats = new ArrayList<>();
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.immigration_index = 0;
        this.tau = tau;

        //here there's been a slight modification so that the deterioration rate now depends on the ratio
        //with the max_growth rate, like in the biofilm_threshold_theory code
        this.deterioration_rate = deterioration_ratio*max_gRate;
        this.biofilm_threshold = biofilm_threshold;

        microhabitats.add(new Microhabitat(K, calc_C_i(0, c_max, alpha, delta_z), scale, sigma, biofilm_threshold));

        microhabitats.get(0).setSurface();
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }

    private BioSystem(double alpha, double c_max, double scale, double sigma, double tau_variable){
        //constructor used to investigate the effects of varying tau step size
        this.alpha = alpha;
        this.c_max = c_max;
        this.scale = scale;
        this.sigma = sigma;
        this.microhabitats = new ArrayList<>();
        this.time_elapsed = 0.;
        this.exit_time = 0.;
        this.immigration_index = 0;
        this.tau = tau_variable;

        this.biofilm_threshold = 0.6;
        this.deterioration_rate = 0.002;

        microhabitats.add(new Microhabitat(K, calc_C_i(0, this.c_max, this.alpha, this.delta_z), scale, sigma, this.biofilm_threshold));
        microhabitats.get(0).setSurface();
        microhabitats.get(0).addARandomBacterium_x_N(5);
    }


    private int getDetachments_counter(){return detachments_counter;}
    private int getDeaths_counter(){return deaths_counter;}
    private int getReplications_counter(){return replications_counter;}
    private int getImmigrations_counter(){return immigrations_counter;}
    private int getTau_halves_counter(){return tau_halves_counter;}

    private double getTimeElapsed(){return time_elapsed;}
    private double getExit_time(){return exit_time;}
    private int getSystemSize(){return microhabitats.size();}


    private int getTotalN(){
        int runningTotal = 0;
        for(Microhabitat m : microhabitats) {
            runningTotal += m.getN();
        }
        return runningTotal;
    }

    private int getBiofilmEdge(){
        int edgeIndex = 0;
        for(int i = 0; i < microhabitats.size(); i++){
            if(microhabitats.get(i).isBiofilm_region()) edgeIndex = i;
        }
        return edgeIndex;
    }

    private int getBiofilmThickness(){
        int thickness = 0;
        for(int i = 0; i < microhabitats.size(); i++){
            if(microhabitats.get(i).isBiofilm_region()) thickness = i+1;
        }
        return thickness;
    }


    private void immigrate(int mh_index, int n_immigrants){
        microhabitats.get(mh_index).addARandomBacterium_x_N(n_immigrants);
    }




    private static double calc_C_i(int i, double c_max, double alpha, double delta_x){
        return c_max*Math.exp(-alpha*i*delta_x);
    }


    public void migrate(int mh_index, int bac_index){

        double migrating_bac = microhabitats.get(mh_index).getPopulation().get(bac_index);
        microhabitats.get(mh_index).removeABacterium(bac_index);

        if(microhabitats.get(mh_index).isSurface()){
            microhabitats.get(mh_index+1).addABacterium(migrating_bac);
        }else if(microhabitats.get(mh_index).isImmigration_zone()){
            microhabitats.get(mh_index-1).addABacterium(migrating_bac);
        }else{
            if(rand.nextBoolean()){
                microhabitats.get(mh_index+1).addABacterium(migrating_bac);
            }else{
                microhabitats.get(mh_index-1).addABacterium(migrating_bac);
            }
        }
    }


    private void updateBiofilmSize(){
        //once the edge microhabitat is sufficiently populated, this adds another microhabitat onto the system list
        //which is then used as the immigration zone

        if(microhabitats.get(immigration_index).atBiofilmThreshold()){

            microhabitats.get(immigration_index).setBiofilm_region();
            microhabitats.get(immigration_index).setImmigration_zone(false);

            int i = microhabitats.size();
            microhabitats.add(new Microhabitat(K, BioSystem.calc_C_i(i, c_max, alpha, delta_z), scale, sigma, biofilm_threshold));
            immigration_index = i;
            microhabitats.get(immigration_index).setImmigration_zone(true);
        }

        //this stops sims going onn unnecessarily too long. if the biofilm reaches the thickness limit then we record the
        //time this happened at and move on
        if(getSystemSize()==thickness_limit){
            exit_time = time_elapsed;
            time_elapsed = 9e9; //this way the time elapsed is now way above the duration value, so the simulation will stop
        }
    }


    public void performAction(){

        double tau_step = tau;

        int system_size = microhabitats.size();
        int[][] replication_allocations;
        int[][] death_allocations;
        int[][] migration_allocations;
        int[] detachment_allocations;
        int[] original_popsizes;
        int n_immigrants;

        whileloop:
        while(true) {
            poiss_imigration = new PoissonDistribution(immigration_rate*tau_step);
            poiss_deterioration = new PoissonDistribution(deterioration_rate*tau_step);
            poiss_migration = new PoissonDistribution(migration_rate*tau_step);
            poiss_migration_edge = new PoissonDistribution(0.5*migration_rate*tau_step);

            replication_allocations = new int[system_size][];
            death_allocations = new int[system_size][];
            migration_allocations = new int[system_size][];
            original_popsizes = new int[system_size];
            detachment_allocations = new int[microhabitats.get(immigration_index).getN()];

            for(int mh_index = 0; mh_index < system_size; mh_index++) {

                //we iterate through all the bacteria and calculate the events which they'll experience
                int mh_pop = microhabitats.get(mh_index).getN();
                int[] n_replications = new int[mh_pop];
                int[] n_deaths = new int[mh_pop];
                int[] n_migrations = new int[mh_pop];

                for(int bac_index = 0; bac_index < mh_pop; bac_index++) {
                    ///////// REPLICATIONS AND DEATHS ///////////////////
                    double[] g_and_d_rate = microhabitats.get(mh_index).replicationAndDeathRates(bac_index);
                    double g_rate = g_and_d_rate[0], d_rate = Math.abs(g_and_d_rate[1]);

                    if(g_rate > 0.) {
                        PoissonDistribution poiss_replication = new PoissonDistribution(g_rate*tau_step);
                        poiss_replication.reseedRandomGenerator(rand.nextLong());
                        n_replications[bac_index] = poiss_replication.sample();
                    }

                    //d_rate is always > 0 due to inclusion of uniform death rate, so no need for the if statements
                    //seen in earlier versions
                    PoissonDistribution poiss_death = new PoissonDistribution(d_rate*tau_step);
                    poiss_death.reseedRandomGenerator(rand.nextLong());
                    n_deaths[bac_index] = poiss_death.sample();

                    //Bacteria cannot replicate and die in the same timestep
                    if(n_deaths[bac_index] > 0 && n_replications[bac_index] > 0){
                        tau_halves_counter++;
                        tau_step /= 2;
                        continue whileloop;
                    }

                    //bacteria can't die twice, so need to handle this
                    if(n_deaths[bac_index] > 1) {
                        tau_halves_counter++;
                        tau_step /= 2;
                        continue whileloop;
                    }


                    ///////// MIGRATIONS AND DETACHMENTS //////////////////////
                    //only non-dead bacteria can migrate or detach
                    if(n_deaths[bac_index] == 0) {

                        //firstly work out the migrations
                        //do edge cases and bulk, then do detachments and set detaching migrations to 0
                        //only do migrations if there's multiple microhabs
                        if(immigration_index > 0) {
                            if(mh_index == 0 || mh_index == immigration_index) {
                                n_migrations[bac_index] = poiss_migration_edge.sample();
                            } else {
                                n_migrations[bac_index] = poiss_migration.sample();
                            }
                            //check for double events
                            if(n_migrations[bac_index] > 1) {
                                //tau_halves_counter++;
                                tau_step /= 2.;
                                continue whileloop;
                            }
                        }

                        //Now do detachments
                        //detaching bacteria can't migrate
                        if(mh_index == immigration_index){
                            detachment_allocations[bac_index] = poiss_deterioration.sample();
                            //check for double events
                            if(detachment_allocations[bac_index] > 1) {
                                //tau_halves_counter++;
                                tau_step /= 2.;
                                continue whileloop;
                            }
                            //bacteria can only migrate if it's not detaching
                            if(detachment_allocations[bac_index] > 0) {
                                n_migrations[bac_index] = 0;
                            }

                        }

                    }
                    //////////////////////////////////////////////////////
                }
                replication_allocations[mh_index] = n_replications;
                death_allocations[mh_index] = n_deaths;
                migration_allocations[mh_index] = n_migrations;
                original_popsizes[mh_index] = microhabitats.get(mh_index).getN();
            }
            n_immigrants = poiss_imigration.sample();
            break whileloop;
        }


        //now we carry out the actions
        for(int mh_index = 0; mh_index < system_size; mh_index++){
            //iterate backwards over the bacteria so we can remove them without getting index errors
            for(int bac_index = original_popsizes[mh_index]-1; bac_index >= 0; bac_index--){

                if(death_allocations[mh_index][bac_index]!= 0) {
                    microhabitats.get(mh_index).removeABacterium(bac_index);
                    deaths_counter++;
                }

                else{
                    microhabitats.get(mh_index).replicateABacterium_x_N(bac_index, replication_allocations[mh_index][bac_index]);
                    replications_counter += replication_allocations[mh_index][bac_index];

                    if(system_size > 1){
                        if(migration_allocations[mh_index][bac_index] != 0) migrate(mh_index, bac_index);
                    }

                    if(mh_index == immigration_index){
                        if(detachment_allocations[bac_index] != 0) {
                            microhabitats.get(mh_index).removeABacterium(bac_index);
                            detachments_counter++;
                        }
                    }
                }
            }
        }

        immigrate(immigration_index, n_immigrants);
        immigrations_counter += n_immigrants;
        updateBiofilmSize();
        //update the time elapsed in the system by the value of tau used in the final events
        time_elapsed += tau_step;

    }


    public static void varyingDeteriorationAndThreshold(int n_reps, double tau_val, int K_startIndex, int K_endIndex){
        long startTime = System.currentTimeMillis();
        //this method varies the deterioration rate and the threshold biofilm density, returns the thickness reached and the event counters
        //this has been modified so that we now vary the deterioration ratio like in the biofilm_threshold_theory notes, rather than
        //the actual rate itself.
        //int n_reps = 15; //the number of times each simulation is repeated for
        //n_reps has been replaced with an argument, which will be the number of processors available to it when I submit this on the qsub system.
        //because this method takes so long to run, we'll split it up into chunks and save the results for each parameter pair seperately
        //the K_start/endIndex arguments tells us which chunks to run, allowing us to run several chunks in parallel and also avoiding the
        //qsub time limit
        int n_measurements = 20; //the number of different values used for deterioration and rho

        double K_min = 0.4, K_max = 1.; //population density threshold for biofilm formation
        double K_increment = (K_max - K_min)/(double)n_measurements;
        double detRatio_min = 0.1, detRatio_max = 0.9;
        double detRatio_increment = (detRatio_max - detRatio_min)/(double)n_measurements;
        double duration = 240.; //10 days


        String directoryName = "/Disk/ds-sopa-personal/s1212500/multispecies-sims/ms_diags_results";
        String filename = String.format("varying-r_det-(%.5f-%.5f)-N_thresh-(%.4f-%.4f)-tau=%.3f", detRatio_min, detRatio_max, K_min, K_max, tau_val);
        String[] headers = new String[]{"tau", "sim_time", "sim_time_stDev","exit_time", "exit_time_stDev", "N*", "det_rate_ratio", "thickness", "thick_stDev", "n_deaths", "n_detachments", "n_immigrations", "n_replications", "n_tau_halves"};
        ArrayList<Databox> Databoxes = new ArrayList<>(); //this isn't really necessary with the new filewriting system

        for(double thresh_K = K_min+(K_startIndex*K_increment); thresh_K <= K_min+(K_endIndex*K_increment); thresh_K+=K_increment){
            for(double det_ratio = detRatio_min; det_ratio <= detRatio_max; det_ratio+=detRatio_increment){
                Databox db = BioSystem.varyingDeteriorationAndThreshold_subroutine(n_reps, duration, thresh_K, det_ratio, tau_val);
                Databoxes.add(db); //not really necessary with the new filewriting system

                //here we can now write the results of each parameter pair to a file.  This will allow our progress to be saved during the simulation
                //seeing as it will take over a week in its current state.
                String solo_filename = String.format("ms_diags-N^-%.3f_rDet-%.5f", thresh_K, det_ratio);
                Toolbox.writeAverageDataboxToFile(directoryName, solo_filename, headers, db);
            }
        }


        Toolbox.writeDataboxArraylistToFile(directoryName, filename, headers, Databoxes);


        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }



    public static Databox varyingDeteriorationAndThreshold_subroutine(int n_reps, double duration, double thresh_K, double det_r_ratio, double tau_val){
        //run several reps of the same parameter set in parallel, then average the results at the end
        Databox[] databoxes = new Databox[n_reps];

        IntStream.range(0, n_reps).parallel().forEach(i -> databoxes[i] = BioSystem.varyingDeteriorationAndThreshold_subsubroutine(i, duration, thresh_K, det_r_ratio, tau_val));

        return Databox.averagedMeasurementsAndStDev(databoxes);
    }




    public static Databox varyingDeteriorationAndThreshold_subsubroutine(int i, double duration, double thresh_K, double det_r, double tau_val){
        int nMeasurements = 50;
        double interval = duration/nMeasurements;
        boolean alreadyRecorded = false;

        BioSystem bs = new BioSystem(det_r, thresh_K, tau_val);
        double start_time = System.currentTimeMillis();

        while(bs.time_elapsed <= (duration+0.001*interval)){
            //got rid of the alreadyRecorded stuff as it doesn't really matter here
            if((bs.getTimeElapsed()%interval >= 0. && bs.getTimeElapsed()%interval <= 0.02*interval)){

                System.out.println("rep : "+i+"\ttau: "+bs.tau+"\tK*: "+bs.biofilm_threshold+"\td_rate: "+bs.deterioration_rate+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+bs.getTotalN()+"\tbf_edge: "+bs.getBiofilmEdge()+"\tsystem size: "+bs.getSystemSize()+"\tc_max: "+bs.c_max);
                //alreadyRecorded = true;
            }

            //if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }

        double finish_time = System.currentTimeMillis();
        double simulation_time = finish_time - start_time;
        double[] counters = new double[]{bs.deaths_counter, bs.detachments_counter, bs.immigrations_counter, bs.replications_counter, bs.tau_halves_counter};

        return new Databox(bs.tau, simulation_time, bs.exit_time, bs.biofilm_threshold, bs.deterioration_rate,  bs.getBiofilmThickness(), counters);
    }






    public static void varyingTauStep(double scale, double sigma){
        //this method was used to investigate the effect of varying the initial tau step
        //you eventually reach a point where the increase in time step is balanced by the number of times you need to
        //half it to avoid double events.
        long startTime = System.currentTimeMillis();
        //this method varies the deterioration rate and the threshold biofilm density, returns the thickness reached and the event counters
        int n_reps = 20; //the number of times each simulation is repeated for
        int n_measurements = 64; //the number of values used taken for tau
        //int n_reps = 4; //the number of times each simulation is repeated for
        //int n_measurements = 3; //the number of values used taken for tau

        double tau_min = 0.01, tau_max = 1.2;
        double tau_increment = (tau_max - tau_min)/(double)n_measurements;
        double duration = 1000.; //1000 hours
        String filename = String.format("varying_tauStep-(%.4f-%.4f)-c=10.0", tau_min, tau_max);
        String[] headers = new String[]{"tau", "sim_time", "sim_time_stDev", "exit_time", "exit_time_stDev", "K*", "det_rate", "thickness", "thick_stDev", "n_deaths", "n_detachments", "n_immigrations", "n_replications", "n_tau_halves"};

        ArrayList<Databox> Databoxes = new ArrayList<>();

        for(double tau = tau_min; tau <= tau_max; tau+=tau_increment){
            Databox db = BioSystem.varyingTauSubroutine(n_reps, duration, scale, sigma, tau);
            Databoxes.add(db);
        }

        Toolbox.writeDataboxArraylistToFile("diagnostics", filename, headers, Databoxes);


        long finishTime = System.currentTimeMillis();
        String diff = Toolbox.millisToShortDHMS(finishTime - startTime);
        System.out.println("results written to file");
        System.out.println("Time taken: "+diff);
    }


    public static Databox varyingTauSubroutine(int n_reps, double duration, double scale, double sigma, double tau){

        Databox[] databoxes = new Databox[n_reps];

        IntStream.range(0, n_reps).parallel().forEach(i -> databoxes[i] = BioSystem.varyingTauSubsubroutine(i, duration, scale, sigma, tau));

        return Databox.averagedMeasurementsAndStDev(databoxes);
    }


    public static Databox varyingTauSubsubroutine(int i, double duration, double scale, double sigma, double tau){

        double c_max = 10.;
        double alpha = 0.01;

        int nMeasurements = 10;
        double interval = duration/nMeasurements;
        boolean alreadyRecorded = false;

        BioSystem bs = new BioSystem(alpha, c_max, scale, sigma, tau);

        double start_time = System.currentTimeMillis();

        while(bs.time_elapsed <= (duration+0.001*interval)){
            if((bs.getTimeElapsed()%interval <= 0.02*interval) && !alreadyRecorded){

                int total_N = bs.getTotalN();
                System.out.println("rep : "+i+"\ttau: "+bs.tau+"\ttau_1/2s: "+bs.tau_halves_counter+"\tK*: "+bs.biofilm_threshold+"\td_rate: "+bs.deterioration_rate+"\tt: "+bs.getTimeElapsed()+"\tpop size: "+total_N+"\tbf_edge: "+bs.getBiofilmEdge()+"\tsystem size: "+bs.getSystemSize()+"\tc_max: "+bs.c_max);
                alreadyRecorded = true;
            }

            if(bs.getTimeElapsed()%interval >= 0.1*interval) alreadyRecorded = false;

            bs.performAction();
        }

        double finish_time = System.currentTimeMillis();
        double simulation_time = finish_time - start_time;

        double[] counters = new double[]{bs.deaths_counter, bs.detachments_counter, bs.immigrations_counter, bs.replications_counter, bs.tau_halves_counter};

        return new Databox(bs.tau, simulation_time, bs.exit_time, bs.biofilm_threshold, bs.deterioration_rate,  bs.getBiofilmThickness(), counters);
    }







}

