This .zip file contains the MATLAB replication codes for Eichenbaum, 
Rebelo and Trabandt (2020), 'The Macroeconomics of Epidemics'.

All codes are written and were run with MATLAB 2016b.

The folder <baseline> contains the code for the basic SIR-macro model. 
Inspect and run sir_macro.m. All other routines in the folder are called automatically.

All other folders in the .zip file contain the same code as in the folder <baseline>, 
yet with different parameter settings. For example, the folder <vaccine> contains the model code 
in which the vaccine parameter deltav is changed from 0 to 1/52.

Likewise, the folder <medprep> contains the code for the medical preparedness model in which 
kappa is set from 0 to 0.9.

The folder <treatment> contains the model in which the treatment parameter deltac is set 
from 0 to 1/52.

The folder <vaccine_treatment_medprep> contains the code for the model named benchmark SIR-Macro model 
in the manuscript. It is the model with deltav=deltac=1/52 and kappa=0.9.

The folders <vaccine_treatment_medprep_optimal_earlyexit_after12w>, <vaccine_treatment_medprep_optimal_earlyexit_after44w> and
<vaccine_treatment_medprep_optimal_lateentry_w33> contain the codes for the early exit and late entry 
scenarios in the benchmark SIR-macro model described in the manuscript. 

The folders ending with _optimal contain the versions of the models in which optimal 
simple containment policy is computed (planner chooses containment rate mu). For example, do_opt_policy=1 is set in sir_macro.m 
in folder <baseline_optimal>. When running sir_macro.m in that folder, optimal policy is computed.
Beware that computations of optimal policy can take some time. We used parallel 
computing with 8 CPU cores on a 3.5 GHZ machine. It took about 10-15 minutes to calculate
optimal policy. On single CPU machines this can take much longer.

The folders ending with _optimal_command_containment compute optimal command containment 
when the planner chooses identical allocations for consumption and hours across health types.
See the paper for the details.

The folder ending with _smart_containment computes the planner solution when the planner chooses
optimal hours and consumption specific to health types (smart containment).

Finally, the files in this .zip file that are named fig1_SIRmacro_SIR.m, fig2_..., fig3_... etc
are the MATLAB .m files that reproduce all the figures provided in the paper. 

Martin Eichenbaum, Sergio Rebelo and Mathias Trabandt (mathias.trabandt@gmail.com)

February, 2021



