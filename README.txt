This program is used for the diagnostic runs of the multispecies simulations in order to find values for the threshold density and deterioration rate parameters.

The performAction() method has now been updated to stop bacteria dying and replicating in the same timestep.

The deterioration rate is now a function of the max gRate, like in the biofilm_threshold_theory notes.

Also changed some of the parameters as discussed.

delta_z: 5 -> 1      (thinner microhabitats)
migration_rate: 0.2 -> 1 (thinner microhabitats)
immigration_ratio: 0.8 -> 80 (new approach to attachment theory)
K: 120 -> 2000 (new microhabitat volume 1mm x 1mm x 1 micron)