// -------------------------------------------------------------------------- //
//                                                                            //
//       Bayesian age-structured model for Prince William Sound herring       //
//                                                                            //
//                               VERSION 1.0                                  //
//                                Jan  2020                                   //
//                                                                            //
//                                 AUTHORS                                    //
//                               John Trochta                                 //                                 
//                              johnt23@uw.edu                                //
//                             Trevor A. Branch                               //
//                              tbranch@uw.edu                                //
//                                                                            //
//                Built on code developed by Melissa Muradian                 //
//           Adapted from Excel-based model by Steven Moffitt (ADF&G)         //
//                                                                            //
// -------------------------------------------------------------------------- //
//                                                                            //
// Program file:  PWS_ASA.tpl                                                 //
// Input Data files                                                           //
//        PWS_ASA.dat:             Model input (surveys, catches, etc.)       //
//        PWS_ASA(phases):         Parameter phases                           //
//        PWS_ASA(ESS).ctl:        Effective sample sizes iteratively         //
//                                 estimated (in R) and used in obj function  //
//        pws_asa.PIN:             Included to test different starting values //
//                                                                            //
// Output files:                                                              //
//        PWS_ASA.rep: Customizable file in the REPORT_SECTION, useful for    //
//                     trouble-shooting                                       //
//        PWS_ASA.std: Default ADMB report file with Hessian derived SE's     //
//        PWS_ASA.par: ADMB parameter estimates                               //
//                                                                            //
// See runADMBmodel.html for explanation on how to run model                  //
//                                                                            //
// -------------------------------------------------------------------------- //


GLOBALS_SECTION
  #include <admodel.h>
  #include <string.h>
  #include <time.h>

  // Following adapted from thread on ADMB users Google Group
  // https://groups.google.com/forum/#!topic/admb-users/WcrSmZc_igw
  dvector rmultinom(const int& seed, const int& size,const dvector& prob)
  {  //Returns a multinomial sample, of size n, based on sampling probabilities p.
  //p is normalized internally, based on the same method employed in R
  random_number_generator rng(seed);
  int i,n,lb,ub;
  float p;
  lb=prob.indexmin(); ub=prob.indexmax();
  dvector freq(lb,ub); freq.initialize();
  dvector P=prob;
  P/=sum(P);
  dvector bisamp(1,size); bisamp.fill_randbi(P[lb],rng);
  freq[lb]=sum(bisamp);
  for(i=lb+1;i<=ub;i++)
  {
  n=size-sum(freq);
  p=P[i]/(1.-sum(P(lb,i-1)));
  //Corrected version

  //cout<<ub-i<<endl;
  dvector bisamp(1,n); bisamp.fill_randbi(p,rng);
  freq[i]=sum(bisamp);
  if(sum(freq)==size) break;
  }
  return (freq);
  }


DATA_SECTION

int DD_Mat;
!! DD_Mat=0;

// |---------------------------------------------------------------------------|
// | CHECK FOR OPTIONAL COMMAND LINE ARGUMENTS & SET FLAGS
// |---------------------------------------------------------------------------|
// | b_simulation_flag  -> flag for running in simulation mode
// | rseed      -> random number seed for simulation 
 int b_simulation_flag;
 int rseed;
 int no_estimation;
 int pin_write;
 LOCAL_CALCS
    int on = 0;
    rseed  = 0;
    no_estimation = 0;
    b_simulation_flag = 0;
    if (ad_comm::argc > 1){
      int on = 0;
      if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-sim")) > -1 ){
        b_simulation_flag = 1;
        //from Merrill: idea for turning on global_parfile flag to look at pin file to fix parameter values for simulation (would need to declare global_parfile and set to zero when declaring b_simulation_flag)
        rseed = atoi(ad_comm::argv[on+1]);
      }

      if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-noest")) > -1 ){
        no_estimation = 1;
      }

      if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-pinwrite")) > -1 ){
        pin_write = 1;
      }
    }

 END_CALCS

// |---------------------------------------------------------------------------|
// | Read in data inputs
// |---------------------------------------------------------------------------|
  //Number of years - nyr
    init_int nrow
    int nyr
    !! nyr=nrow;

  // Number of years to be fit in the model
    init_int nyr_tobefit
 
  //Number of age classes - 7
    init_int ncol
    int nage
    !! nage=ncol;

  //Weight-at-age
    init_matrix w_a_a(1,nyr,1,nage)
  
  //Fecundity-at-age
    init_matrix fecun(1,nyr,1,nage)
  
  //Pound Catch
    init_matrix pc(1,nyr,1,nage)
  
  //Proportion of pound fish killed
    init_number pk
  
  //Food/Bait Catch
    init_matrix fbc(1,nyr,1,nage)
  
  //Gillnet catch
    init_matrix gc(1,nyr,1,nage)
  
  //Total seine yield
    init_vector sc(1,nyr)
  
  //% Female Spawners
    init_vector f_sp(1,nyr)
  
  //MDM
    init_vector mdm(1,nyr)
    
  //Egg Deposition
    init_vector egg(1,nyr)
    
  //Standard errors derived from confidence intervals given for Egg Deposition ln(s.e.)
    init_vector cv_egg(1,nyr)       
    
  //ADFG Hydroacoustic Survey data - is a combination of ADFG and PWSSC estimates until 1994, from '95 onward are only ADFG survey data
    init_number hydADFG_start //read in year from the data file as first year of the survey: 1995
    init_vector hydADFG(1,nyr)
    
  //PWSSC Hydroacoustic Survey data 
    init_number hydPWSSC_start //read in year  from the data file as first year of the survey: 1993
    init_vector hydPWSSC(1,nyr)
    
  //Standard errors derived from confidence intervals given for PWSSC hydroacoustic biomass ln(s.e.)
    init_vector cv_hydPWSSC(1,nyr)

  //Seine age distribution
    init_matrix seine(1,nyr,1,nage)

  //Spawning age composition
    init_matrix spac(1,nyr,1,nage)

  // Aerial juvenile survey (incorporated 12/2019)
  init_vector juv_ind(1,nyr)

  //Maturity model - selects which one to use in the ASA
  // 1 is the current model, 2 is model incorporating spawner survey selectivity and fit to maturity indices
    init_int mat_mod_type

  // Seroprevalence: nyr x 2*nage matrix with each column representing # of positive then negative samples in each age (e.g. # positive age 0, # negative age 0, etc.)
    int n_sero
    !! n_sero=2*nage;
    init_matrix sero_obs(1,nyr,1,n_sero)  
    init_number sero_start    //start year of observations for fitting
    init_number vhsv_start_est    //start year for estimating infection (can be earlier than sero_start)
    init_int recov_prob_type

// |---------------------------------------------------------------------------|
// | Read in parameter phases
// |---------------------------------------------------------------------------|
  !! ad_comm::change_datafile_name("PWS_ASA(phases).ctl");   // has differing phases up to 5

    init_int ESS_est   
    init_int ph_Z_0_8
    init_int ph_Z_9
    init_int ph_Z_0_8offset
    init_int ph_Z_9offset
    init_int ph_matur_age3_per1
    init_int ph_matur_age4_per1
    init_int ph_matur_age3_per2
    init_int ph_matur_age4_per2
    init_int ph_alpha_v
    init_int ph_beta_v
    init_int ph_age0devs
    init_int ph_init_pop
    init_int ph_eggAdd
    init_int ph_mdm
    init_int ph_mdmAdd
    init_int ph_hyd1
    init_int ph_hydAdd1
    init_int ph_hyd2
    init_int ph_hydAdd2
    init_int ph_age3_4mort_93
    init_int ph_age5_8mort_93
    
    init_int ph_meanage0
    init_int ph_meanage0_offset
    init_int ph_sigmaage0
    init_int ph_betaage0  
    init_int ph_betamortality       
    init_int ph_mortdevs

    init_int ph_age0_offset
    init_int ph_mortality_offset

    init_int ph_survey_vul_alpha
    init_int ph_survey_vul_beta

    init_int ph_sigma_mortdevs
    init_int ph_sigma_age0covar
    init_int ph_sigma_morcovar

    // Aerial juvenile survey (incorporated 12/2019)
    init_int ph_juv_ind

    // Seroprevalence parameters
    init_int ph_vhs_pars
    init_int ph_vul_symp
    init_int ph_vul_sero

// |---------------------------------------------------------------------------|
// | Read in effective sample sizes for age-composition (calculated reiteratively before running MCMC)
// |---------------------------------------------------------------------------|
  !! if (ESS_est == -1) {
  !! ad_comm::change_datafile_name("PWS_ASA(ESS).ctl");}
  !! else if (ESS_est == 1) {
  !! ad_comm::change_datafile_name("PWS_ASA(ESS_estimate).ctl");}

    // Seine Effective Sample Size - these are filled in by R at the end of the PWS_ASA(Data).ctl
      init_vector ESS_Se(1,nyr_tobefit)
    // Spawing Effective Sample Size
      init_vector ESS_Sp(1,nyr_tobefit)
    // Antibody sample size
      init_vector ESS_Antibody(1,nyr_tobefit)

// |---------------------------------------------------------------------------|
// | Read in the recruitment and natural mortality deviate information
// |---------------------------------------------------------------------------|
  !! ad_comm::change_datafile_name("PWS_ASA(covariate).ctl");  

    init_int standardize_covariates
  	init_int n_age0_covs
    init_ivector R_fixed_or_est(1,n_age0_covs)
    init_ivector age0_turn_on(1,n_age0_covs)
    init_matrix age0_covariates(1,nyr,1,n_age0_covs)
    init_ivector R_change(1,nyr)

    init_int n_mor_covs
    init_ivector M_fixed_or_est(1,n_mor_covs)
    init_ivector mor_season(1,n_mor_covs)
    init_ivector mor_turn_on(1,n_mor_covs)
    init_matrix covariate_effect_byage(1,nage,1,n_mor_covs)
    init_matrix mor_covariates(1,nyr,1,n_mor_covs)
    init_vector nyr_tobefit_winter_covariate(1,n_mor_covs)
    init_ivector M_change(1,nyr)

// |---------------------------------------------------------------------------|
// | Read-in settings for simulations
// |---------------------------------------------------------------------------|
  !! ad_comm::change_datafile_name("PWS_ASA(sim_settings).ctl");   

    // I have included an option for user to input effort to simulate catches - this turns it on/off
    init_int sim_catches
    init_int nyr_resample_period
    init_ivector resample_period(1,nyr_resample_period)
    init_vector exploitation_history(1,nyr_tobefit)
    init_int age0_dev_option
    init_int data_avg_option

    // Switches to 1 within simulations function if sim_catches=1
    int turn_on_effort
    !! turn_on_effort=0;

    matrix gc_V(1,nyr_tobefit,1,nage)
    matrix pc_V(1,nyr_tobefit,1,nage)
    matrix fbc_V(1,nyr_tobefit,1,nage)
    vector sc_F(1,nyr_tobefit)

    //!!cout << "Done Data" << endl;

// |---------------------------------------------------------------------------|
// | Prelim value assignments & writing the PIN file if specified
// |---------------------------------------------------------------------------|
   
   // Lower and upper bounds to maturity ogive proportions
   number LB_Mat3_1
   number LB_Mat3_2
   number LB_Mat4_1
   number LB_Mat4_2
   number UB_Mat3_1
   number UB_Mat3_2
   number UB_Mat4_1
   number UB_Mat4_2

   // This sets up the estimable recruit and mortality covariate parameters
   // If mortality covariates are not included, counter is set to one so program does not crash (used for specifying size of estimable vector)
   // If mortality covariates are included, beta and/or deviate parameters on mortality are turned ON
   int mor_cov_counter
   int M_cov_model
 
   // M_cov_model default is 1, whether or not covariates are included in the model
   // M_cov_model changes to 2 if even just one of the covariates is being modeled as
   // an index that is a prior in the model
   !! M_cov_model=1;
   !! for(int i=1; i<=n_mor_covs; i++){
   !!   if(M_fixed_or_est(i)*mor_turn_on(i)==2){
   !!     M_cov_model=2;
   !!   }
   !! }
   !! mor_cov_counter=sum(mor_turn_on);
   !! if(mor_cov_counter==0){
   !!   mor_cov_counter=1;
   !! }else{
   !!   ph_betamortality=2;
   !!   if(M_cov_model==2){
   !!     ph_mortdevs=2;
   !!   }
   !! }

   int rec_cov_counter
   int rec_cov_counter_age0devs
   int R_cov_model
   number sigma_age0devs_PIN

   // R_cov_model default is 1, whether or not covariates are included in the model
   // R_cov_model changes to 2 if even just one of the covariates is being modeled as
   // an index that is a prior in the model
   !! R_cov_model=1;
   !! for(int i=1; i<=n_age0_covs; i++){
   !!   if(R_fixed_or_est(i)*age0_turn_on(i)==2){
   !!     R_cov_model=2;
   !!   }
   !! }
   !! rec_cov_counter=sum(age0_turn_on);
   !! if(rec_cov_counter==0){
   !!   rec_cov_counter=1;
   !! }else{
   !!   ph_betaage0=2;
   !! }
   !! sigma_age0devs_PIN=0;
   !! if(R_cov_model==1){
   !!   rec_cov_counter_age0devs=1;
   !! }else if(R_cov_model==2){
   !!   rec_cov_counter_age0devs=rec_cov_counter;
   !! }

   // Other PIN parameters
   number Mat3_1_PIN
   number Mat3_2_PIN
   number Mat4_1_PIN
   number Mat4_2_PIN
   number sigma_mortdevs_PIN

   vector sigma_age0covar_PIN(1,rec_cov_counter_age0devs) // Weight of weighted Sum of Squares fit to recruitment indices
   vector sigma_morcovar_PIN(1,mor_cov_counter) // Weight of weighted SS fit to mortality indices
   
   vector beta_age0_PIN(1,rec_cov_counter)
   vector beta_mortality_PIN(1,mor_cov_counter)
   matrix annual_age0devs_PIN(1,rec_cov_counter_age0devs,1,nyr_tobefit-3)
   matrix annual_mortdevs_PIN(1,mor_cov_counter,1,nyr_tobefit)
   vector beta_age0_offset_PIN(1,rec_cov_counter)
   vector beta_mortality_offset_PIN(1,mor_cov_counter)

   vector beta_mortality_ind(1,mor_cov_counter)
   vector beta_recruit_ind(1,rec_cov_counter)

   // VHSV seroprevalence related parameters (11/2020)
   vector annual_inf_PIN(1,nyr_tobefit);
   vector recov_prob_PIN(1,nyr_tobefit);
   number infec_vul_a50_PIN;
   number infec_vul_a95_PIN;
   number sero_samp_vul_a50_PIN;
   number sero_samp_vul_a95_PIN;
   

 LOCAL_CALCS
   int j=1;
   beta_age0_PIN=0;
   annual_age0devs_PIN=0;
   beta_age0_offset_PIN=0;
   for(int i=1; i<=n_age0_covs; i++){
      if(age0_turn_on(i)==1){
        beta_recruit_ind(j)=i;
        beta_age0_PIN(j)=0.1;
        j+=1;
      }
   }
   
   beta_mortality_PIN=0;
   annual_mortdevs_PIN=0;
   beta_mortality_offset_PIN=0;
   j=1;
   for(int i=1; i<=n_mor_covs; i++){
      if(mor_turn_on(i)==1){
        beta_mortality_ind(j)=i;
        beta_mortality_PIN(j)=0.1;
        j+=1;

      }

   }

   sigma_mortdevs_PIN=1;

   sigma_age0covar_PIN=0.5;
   sigma_morcovar_PIN=0.5;

   if(mat_mod_type==1){

     Mat3_1_PIN=0.5;
     Mat4_1_PIN=0.8;
     Mat3_2_PIN=0.5036630;
     Mat4_2_PIN=0.9;

     // Mat3_1_PIN=0.6562465;
     // Mat4_1_PIN=0.8137432;
     // Mat3_2_PIN=0.5036630;
     // Mat4_2_PIN=0.9;
   }else{
     Mat3_1_PIN=0.6562465;
     Mat4_1_PIN=0.8137432;
     Mat3_2_PIN=0.5036630;
     Mat4_2_PIN=0.3818182;

     ph_matur_age3_per1=-3;
     ph_matur_age4_per1=-3;
     ph_matur_age3_per2=-3;
     ph_matur_age4_per2=-3;

     ph_survey_vul_alpha=4;
     ph_survey_vul_beta=4;
   }
   
   // VHSV seroprevalence parameter initial values (11/2020)
   annual_inf_PIN = 0;
   recov_prob_PIN = 0.7;
   infec_vul_a50_PIN = 1;
   infec_vul_a95_PIN = 2;
   sero_samp_vul_a50_PIN = 2;
   sero_samp_vul_a95_PIN = 3;

   if(pin_write){
      ofstream write_pin("PWS_ASA.PIN",ios::trunc);
         
      write_pin << 0 << endl;
      write_pin << 0 << endl;
      write_pin << 0.25 << endl;
      write_pin << 0.827 << endl;
      write_pin << 0 << endl;
      write_pin << 0 << endl;
      write_pin << Mat3_1_PIN << endl;
      write_pin << Mat4_1_PIN << endl;
      write_pin << Mat3_2_PIN << endl;
      write_pin << Mat4_2_PIN << endl;
      write_pin << 3.97004158704 << endl;
      write_pin << 2.37479158941 << endl;
      write_pin << 4 << endl;
      write_pin << 2.4 << endl;
      write_pin << "6.35103457941 5.66604017924 5.92140625419 6.7410980789 4.73590293732"<< endl;
      write_pin << 0.4 << endl;
      write_pin << 5.98 << endl;
      write_pin << 0.305 << endl;
      write_pin << -0.466 << endl;
      write_pin << 0.249 << endl;
      write_pin << -0.387 << endl;
      write_pin << 0.305 << endl;

      // Aerial juvenile survey (incorporated 12/2019)
      write_pin << -1 << endl;
      write_pin << 0 << endl;

      write_pin << annual_age0devs_PIN << endl;
      write_pin << 6.20555781195 << endl;
      write_pin << 0 << endl;
      write_pin << sigma_age0devs_PIN << endl;

      write_pin << beta_age0_PIN << endl;
      write_pin << beta_mortality_PIN << endl;
      write_pin << annual_mortdevs_PIN << endl;
      write_pin << sigma_mortdevs_PIN << endl;
      write_pin << beta_age0_offset_PIN << endl;
      write_pin << beta_mortality_offset_PIN << endl;

      write_pin << sigma_age0covar_PIN << endl;
      write_pin << sigma_morcovar_PIN << endl;

      write_pin << annual_inf_PIN << endl;
      write_pin << recov_prob_PIN << endl;
      write_pin << infec_vul_a50_PIN << endl;
      write_pin << infec_vul_a95_PIN << endl;
      write_pin << sero_samp_vul_a50_PIN << endl;
      write_pin << sero_samp_vul_a95_PIN << endl;
   }

   // LB_Mat3_1=0; 
   // LB_Mat3_2=0;
   // LB_Mat4_1=0.2;
   // LB_Mat4_2=0.2;
   // UB_Mat3_1=0.95;
   // UB_Mat3_2=0.95;
   // UB_Mat4_1=1;
   // UB_Mat4_2=1;

   // For logistic function form of maturity function
   LB_Mat3_1=2; 
   LB_Mat3_2=2;
   LB_Mat4_1=3;
   LB_Mat4_2=3;
   UB_Mat3_1=4;
   UB_Mat3_2=4;
   UB_Mat4_1=6.5;
   UB_Mat4_2=6.5;
   
   if(DD_Mat==1){
    LB_Mat3_1=0; 
    LB_Mat3_2=-12;
    LB_Mat4_1=-4;
    LB_Mat4_2=0.2;
       
    UB_Mat3_1=3;
   	UB_Mat3_2=-7;
   	UB_Mat4_1=4;
   	UB_Mat4_2=1;
   }
 END_CALCS


PARAMETER_SECTION

// |---------------------------------------------------------------------------|
// | NATURAL MORTALITY PARAMETERS
// |---------------------------------------------------------------------------|
  //Age-related instantaneous mortality
  init_bounded_number VHSV_age3_4_mort_93(-1,5,ph_age3_4mort_93)
  init_bounded_number ICH_age5_8_mort_93(-1,5,ph_age5_8mort_93)

  //Estimate baseline adult mortality in log-space - very non-informative prior
  init_bounded_number Z_0_8(0.05,2.3,ph_Z_0_8)        
  
  //implements a constraint from ADF&G model: .25 <= S_9+ <= .95*S_5+, which 
  //than age3-8 natural mortality says that plus group mortality must be larger
  init_bounded_number Z_9(0.30,1.4,ph_Z_9) 

  init_bounded_number Z_0_8offset(-0.5,2,ph_Z_0_8offset)   
  init_bounded_number Z_9offset(-0.25,2,ph_Z_9offset)        
 
// |---------------------------------------------------------------------------|
// | MATURITY PARAMETERS
// |---------------------------------------------------------------------------|
  //Maturity parameters of age 3 and 4 before (per1) and after 1997 (per2)
  init_bounded_number matur_age3_per1(LB_Mat3_1,UB_Mat3_1,ph_matur_age3_per1)
  init_bounded_number matur_age4_per1(LB_Mat4_1,UB_Mat4_1,ph_matur_age4_per1)
  init_bounded_number matur_age3_per2(LB_Mat3_2,UB_Mat3_2,ph_matur_age3_per2)
  init_bounded_number matur_age4_per2(LB_Mat4_2,UB_Mat4_2,ph_matur_age4_per2)

// |---------------------------------------------------------------------------|
// | SELECTIVITY PARAMETERS
// |---------------------------------------------------------------------------|
  //Seine Vulnerability parameters
  init_bounded_number alpha_v(3,5,ph_alpha_v)        
  init_bounded_number beta_v(1,7,ph_beta_v)

  //Survey Vulnerability parameters
  init_bounded_number survey_vul_alpha(3,5,ph_survey_vul_alpha)        
  init_bounded_number survey_vul_beta(1,7,ph_survey_vul_beta)          
  
  //Initial abundance parameters (age-3 for all yrs, 1980 ages 4 and 5+)
  init_bounded_vector loginit_pop(1,5,3,8,ph_init_pop)

// |---------------------------------------------------------------------------|
// | SURVEY SCALAR & CV PARAMETERS
// |---------------------------------------------------------------------------|
  // Egg deposition additional variance term
  init_bounded_number egg_add(0.00001,0.5,ph_eggAdd) // In the likelihoods I use sqrt(), so bound must be
                                                     // positive, since is non-differentiable at 0...

  //Milt coefficient
  init_bounded_number logmdm_c(2.3,7,ph_mdm)         // infinity if mdm_c goes to zero
  init_bounded_number m_add(0.00001,0.9,ph_mdmAdd)

  //Hydroacoustic scalers and additional variance parameters
  init_bounded_number hydADFG_q(-5,5,ph_hyd1) 
  init_bounded_number hydADFG_add(0.00001,0.7,ph_hydAdd1)

  init_bounded_number hydPWSSC_q(-5,5,ph_hyd2) 
  init_bounded_number hydPWSSC_add(0.00001,0.6,ph_hydAdd2)

  // Aerial juvenile survey (incorporated 12/2019)
  init_bounded_number log_q_juv(-5,5,ph_juv_ind)
  init_bounded_number log_overdisp_juv(-5,5,ph_juv_ind)
 
// |---------------------------------------------------------------------------|
// | REC DEVIATES & COVARIATE EFFECTS ON MORTALITY
// |---------------------------------------------------------------------------| 
  // Covariate effects on recruits (age 0) or mortality
  init_bounded_matrix annual_age0devs(1,rec_cov_counter_age0devs,1,nyr_tobefit-3,-10,10,ph_age0devs)
  init_bounded_number log_MeanAge0(2,8,ph_meanage0)
  init_bounded_number Mean_Age0offset(-2,2,ph_meanage0_offset)
  init_bounded_number sigma_age0devs(0.0000001,3,ph_sigmaage0)
  init_bounded_vector beta_age0(1,rec_cov_counter,-20,20,ph_betaage0)
  
  init_bounded_vector beta_mortality(1,mor_cov_counter,-50,50,ph_betamortality)
  init_bounded_matrix annual_mortdevs(1,mor_cov_counter,1,nyr_tobefit,-10,10,ph_mortdevs) // Just estimate the deviates for 1992-1995
  init_bounded_number sigma_mortdevs(0.0000001,3,ph_sigma_mortdevs)

  // Basically sets 2nd time block for different effect estimate on covariates - need to improve....
  init_bounded_vector beta_age0_offset(1,rec_cov_counter,-5,5,ph_age0_offset)
  init_bounded_vector beta_mortality_offset(1,mor_cov_counter,-5,5,ph_mortality_offset)

  init_bounded_vector sigma_age0covar(1,rec_cov_counter_age0devs,0.00001,2,ph_sigma_age0covar) // Weight of weighted Sum of Squares fit to recruitment indices
  init_bounded_vector sigma_morcovar(1,mor_cov_counter,0.00001,2,ph_sigma_morcovar) // Weight of weighted SS fit to mortality indices

// |---------------------------------------------------------------------------|
// | AGE-STRUCTURED DISEASE EFFECTS BASED ON SEROPREVALENCE DATA
// |---------------------------------------------------------------------------|
  init_bounded_vector annual_inf(1,nyr_tobefit,0,1,ph_vhs_pars) 
  init_bounded_vector recov_prob(1,nyr_tobefit,0.2,0.9,ph_vhs_pars)

  init_bounded_number infec_vul_a50(-1,4.5,ph_vul_symp)        
  init_bounded_number infec_vul_a95(0,8,ph_vul_symp) 

  init_bounded_number sero_samp_vul_a50(0,4,ph_vul_sero)        
  init_bounded_number sero_samp_vul_a95(3,6,ph_vul_sero) 

// |---------------------------------------------------------------------------|
// | CALCULATED NUMBERS
// |---------------------------------------------------------------------------|
  number S_0_2 							// Survival rate of ages 0-2
  number S_3_8						  // Survival rate of ages 3-8
  number S_9							  // Survival rate of 9+ age group

  number mdm_c 
  number MDMtemp_2
  number HtempADFG_num

  number not_below_this

 //To count instances of incurring the posfun penalty
  number penCount

  number Se_llk
  number Sp_llk
  number MDMllk
  number EGGllk
  number H_ADFGllk
  number H_PWSSCllk
  number age0_devs_penllk
  number mort_devs_penllk
  number age0_covar_prior
  number mort_covar_prior
  number Z_prior
  number hydADFG_add_prior
  number hydPWSSC_add_prior
  number m_add_prior
  
  // Aerial juvenile survey (incorporated 12/2019)
  number juv_llk

  number temp1MeanLgRec
  number temp2MeanLgRec
  number meanLgRec
  number projected_PFRB //is an sdreport variable

  // Seroprevalence likelihood
  number sero_llk

// |---------------------------------------------------------------------------|
// | CALCULATED VECTORS
// |---------------------------------------------------------------------------|
  vector init_age_0(1,nyr_tobefit)       //Estimates of Age 0's (millions)
  vector forecast_winter_effect(1,nage)
  vector forecast_Sur_winter(1,nage)
  vector age0_effect(1,nyr_tobefit)

  vector Vul(1,nage)                     //Vulnerability of fish to seine commercial fishery
  vector Vul_survey(1,nage)              //Vulnerability of fish to seine gear used for surveys
  vector CC_bot(1,nyr_tobefit)           //Total Available population for Seine catch for each year
  vector N_se(1,nyr_tobefit)             //Seine-Catch, Estimated total catch
  vector SSB(1,nyr_tobefit)              //Pre-Fishery Run Biomass, mt
  vector SB_star(1,nyr_tobefit)          //Management reference point - is SSB (pre-fishery run biomass) in tons
  vector SB(1,nyr_tobefit)               //Total Post-Fishery Spawning Biomass
  vector SpAC_bot(1,nyr_tobefit)         //Total Pre-fishery spawning pop (mil)
  vector SeBiomass(1,nyr_tobefit)        //Sum(Weight-at-age*Seine age comp) is like Seine biomass per year
  vector Early_biomass(1,nyr_tobefit)    //Pre-Fishery Biomass by year

  // Aerial juvenile survey (incorporated 12/2019)
  vector juv_pred(1,nyr_tobefit)

  vector MDM(1,nyr_tobefit)              //Mile-days of milt
  vector EGG(1,nyr_tobefit)              //Egg Deposition - rowsums of EggAC
  vector EggAC_sum(1,nyr_tobefit)        //totol egg dep (both male and female) by year
  vector HYD_ADFG(1,nyr_tobefit)         //ADFG Hydroacoustic Estimates
  vector HYD_PWSSC(1,nyr_tobefit)        //PWSSC Hydroacoustics Estimates
  
  vector MDMtemp_1(1,nyr_tobefit)         //Mile-days of milt
  vector EGGtemp(1,nyr_tobefit)           //Egg deposition
  vector HtempADFG_vec(1,nyr_tobefit)     // ADFG Hydroacoustics
  vector HtempPWSSC_vec(1,nyr_tobefit)    // PWSSC Hydroacoustics
  vector Setemp_2(1,nyr_tobefit)
  vector Sptemp_2(1,nyr_tobefit)
  vector Setemp_3(1,nyr_tobefit)
  vector Sptemp_3(1,nyr_tobefit)
  
 //Analytical Sigmas
  vector Eg_SD(1,nyr_tobefit)
  vector PWSSC_SD(1,nyr_tobefit)
 
 //Likelihood components 
  vector MDMllk_ind(1,nyr_tobefit)
  vector EGGllk_ind(1,nyr_tobefit)
  vector H_ADFGllk_ind(1,nyr_tobefit)
  vector H_PWSSCllk_ind(1,nyr_tobefit)

  // Aerial juvenile survey (incorporated 12/2019)
  vector juv_llk_ind(1,nyr_tobefit)

 // Variables and vectors for calculating projected final year biomass
  vector tempWgt(1,nage)
  vector avgWgt5Yr(1,nage)
  vector projected_N_y_a(1,nage)
  vector projected_Early_Sp_biomass(1,nage)

  vector Mean_Age0(1,nyr_tobefit) 
  vector init_pop(1,5)
  vector agggregate_annual_age0devs(1,nyr_tobefit)
  vector temp_age0(1,rec_cov_counter)

  // Age-specific vulnerability to viral transmission (i.e. proportion mixing with sympatric pop)
  vector Vul_sympat(1,nage)
  vector Vul_sero_survey(1,nage)

  vector inf_inc_sp(1,nyr_tobefit)
  vector seroprev_sp(1,nyr_tobefit)
  vector fatal_sp(1,nyr_tobefit)

  vector inf_inc_age3(1,nyr_tobefit)
  vector seroprev_age3(1,nyr_tobefit)
  vector fatal_age3(1,nyr_tobefit)  

// |---------------------------------------------------------------------------|
// | CALCULATED MATRICES
// |---------------------------------------------------------------------------|
 //Mortality and survival rates
  matrix summer_effect(1,nyr_tobefit,1,nage)
  matrix winter_effect(1,nyr_tobefit,1,nage)
  matrix Sur_summer(1,nyr_tobefit,1,nage)     //Half-year Survival between spring and fall
  matrix Sur_winter(1,nyr_tobefit,1,nage)     //Half-year Survival between fall and spring

 //Age matrices
  matrix Mat(1,nyr_tobefit,1,nage)            //Maturity 
  matrix N_y_a(1,nyr_tobefit,1,nage)          //Pre-Fishery Abundance (millions)
  matrix CC_top(1,nyr_tobefit,1,nage)         //Vulnerability_Seine*N_y_a, avail pop at age
  matrix SeAC(1,nyr_tobefit,1,nage)           //Seine Age Composition
  matrix N_sp(1,nyr_tobefit,1,nage)           //Spawning Population (millions)
  matrix Early_Sp_biomass(1,nyr_tobefit,1,nage)
  matrix SPB(1,nyr_tobefit,1,nage)            //Spawning Biomass at age
  matrix SpAC_top(1,nyr_tobefit,1,nage)       //Pre-Fishery spawning population (millions)
  matrix SpAC(1,nyr_tobefit,1,nage)           //Spawning Age Composition
  matrix bio_a_a(1,nyr_tobefit,1,nage)        //Weight-at-age*Seine age comp, is like Seine biomass-at-age
  matrix EggAC(1,nyr_tobefit,1,nage)          //Egg dep by age, like Egg deposition Age-Comp  
  matrix Early_bio_a_a(1,nyr_tobefit,1,nage)  //N_y_a*w_a_a, Pre-Fishery biomass-at-age

 //Temporary matrices/vectors for calculating likelihood components
  matrix SeACR(1,nyr_tobefit,1,nage)      //Seine age comp
  matrix SpACR(1,nyr_tobefit,1,nage)      //Spawning age comp
  matrix Setemp_1(1,nyr_tobefit,1,nage)
  matrix Sptemp_1(1,nyr_tobefit,1,nage)

  matrix annual_mortdevs_byage(1,nyr_tobefit,1,nage)

  matrix Mat_unobs(1,nyr_tobefit,1,nage)      //Maturity of unobserved population

  // Age-specific matrices related to epidemiological stages
  matrix suscep(1,nyr_tobefit,1,nage)
  matrix immune(1,nyr_tobefit,1,nage)
  matrix Sur_vhsv(1,nyr_tobefit,1,nage)
  matrix sero_pred(1,nyr_tobefit,1,n_sero)

// |---------------------------------------------------------------------------|
// | OBJECTIVE FUNCTION
// |---------------------------------------------------------------------------|
  objective_function_value f_llk
  
  sdreport_number SSB_final_year
  //sdreport_vector init_pop(1,5)         
  


PRELIMINARY_CALCS_SECTION
  if(standardize_covariates){

  // Standardize Age 0 first
    for(int k=1; k<=n_age0_covs; k++){

      double varsums=0.0;
      double varN=0.0;
      double varmean=0.0;
      double varSD=0.0;
      for(int i=1; i<=nyr_tobefit; i++){
        if(age0_covariates(i,k)!=-9){
          varsums+=age0_covariates(i,k);
          varN+=1;
        }
      }
      varmean=varsums/varN;
      for(int i=1; i<=nyr_tobefit; i++){
        if(age0_covariates(i,k)!=-9){
          varSD+=square(age0_covariates(i,k)-varmean);
        }
      }
      varSD=sqrt(varSD/(varN-1));
      for(int i=1; i<=nyr_tobefit; i++){
        if(age0_covariates(i,k)!=-9){
          age0_covariates(i,k)=(age0_covariates(i,k)-varmean)/varSD;
        }
      }
    }

  // Standardize mortality indices 
    for(int k=1; k<=n_mor_covs; k++){
      double varsums=0.0;
      double varN=0.0;
      double varmean=0.0;
      double varSD=0.0;
      for(int i=1; i<=nyr_tobefit; i++){
        if(mor_covariates(i,k)!=-9){
          varsums+=mor_covariates(i,k);
          varN+=1;
        }
      }
      varmean=varsums/varN;
      for(int i=1; i<=nyr_tobefit; i++){
        if(mor_covariates(i,k)!=-9){
          varSD+=square(mor_covariates(i,k)-varmean);
        }
      }
      varSD=sqrt(varSD/(varN-1));
      for(int i=1; i<=nyr_tobefit; i++){
        if(mor_covariates(i,k)!=-9){
          mor_covariates(i,k)=(mor_covariates(i,k)-varmean)/varSD;
        }
      }

      if(nyr_tobefit_winter_covariate(k)!=-9){
        nyr_tobefit_winter_covariate(k)=(nyr_tobefit_winter_covariate(k)-varmean)/varSD;
      }
    }

  }

  if(b_simulation_flag && rseed >= 0){

    cout<<"|--------------------------|"<<endl;
    cout<<"| RUNNING SIMULATION MODEL |"<<endl;
    cout<<"|--------------------------|"<<endl;
    runSimulationModel(rseed);
  
  }
  if(no_estimation){  
    Mean_Age0 = exp(log_MeanAge0+Mean_Age0offset*R_change(1,nyr_tobefit));
    mdm_c = exp(logmdm_c);
    init_pop = exp(loginit_pop);

    calc_naturalmortality();
    // cout << summer_effect_tofit << endl;

    if(DD_Mat==0){
	   calc_maturity();
  	}
    calc_selectivity();
    calc_statevariables();
    calc_surveyvalues();
    calc_nll_components();  

    ofstream deterministic_run("deterministic_run.rep",ios::trunc);
    
    // Aerial juvenile survey (incorporated 12/2019)
    deterministic_run << "# Posterior Probability" << endl;
    deterministic_run << Se_llk +Sp_llk +EGGllk +H_ADFGllk +H_PWSSCllk +MDMllk +age0_devs_penllk +mort_devs_penllk +Z_prior +hydADFG_add_prior +hydPWSSC_add_prior +m_add_prior +juv_llk +sero_llk << endl << endl;

    // Aerial juvenile survey (incorporated 12/2019)
    deterministic_run << "# Likelihood Components & Priors" << endl;
    deterministic_run << "# Se_llk  Sp_llk  EGGllk  H_ADFGllk  H_PWSSCllk  MDMllk  age0_devs_penllk  mort_devs_penllk   Z_prior hydADFG_add_prior hydPWSSC_add_prior  m_add_prior  juv_llk  sero_llk" << endl;
    deterministic_run << Se_llk << "  " << Sp_llk << "  " << EGGllk << "  " << H_ADFGllk << "  " << H_PWSSCllk << "  " << MDMllk << "  " << age0_devs_penllk << "  " << mort_devs_penllk << "  " << Z_prior << " " << hydADFG_add_prior << " " << hydPWSSC_add_prior << " " <<  m_add_prior << " " <<  juv_llk <<  " " <<  sero_llk << endl << endl;

    deterministic_run << "# Pre-fishery Spawning Biomass (metric tons)" << endl;
    deterministic_run << SSB << endl << endl;
    deterministic_run << "# Post-fishery Spawning Biomass (metric tons)" << endl;
    deterministic_run << SB << endl << endl;
    deterministic_run << "# Recruitment (millions of Age 3 herring)" << endl;
    for(int i=1; i<=nyr_tobefit; i++){
      deterministic_run << N_y_a(i,4) << " ";
    }
    deterministic_run << endl;
  }
  // cout << "complete preliminary" << endl;

PROCEDURE_SECTION

// |---------------------------------------------------------------------------|
// | MODEL ROUTINES
// |---------------------------------------------------------------------------|
// | PSUEDOCODE:
// | - initialize likelihood and penalty counts (e.g. model calculates negative biomass or >100% survival)
// | - calculate survival from natural mortality
// | - calculate maturity ogives
// | - calculate fishery selectivity
// | - Calculate state variables:
// |    - calculate spawning stock biomass & age composition
// | - observation models:
// |    - create observed values of egg deposition, milt-days, acoustic, & age-comp survey
// | - calculate objective function value
// | - evaluate additional functions in mceval_phase (biomass forecase & writing result files)
// |---------------------------------------------------------------------------|
  
  Mean_Age0 = exp(log_MeanAge0+Mean_Age0offset*R_change(1,nyr_tobefit));
  mdm_c = exp(logmdm_c);
  init_pop = exp(loginit_pop);
  f_llk=0;
  penCount = 0;

  calc_naturalmortality();

  if(DD_Mat==0){
	calc_maturity();
  }

  calc_selectivity();

  calc_statevariables();

  calc_surveyvalues();
  
  calc_nll_components();

  // Objective function ADMB will minimize
  // Aerial juvenile survey (incorporated 12/2019)
  f_llk += Se_llk +Sp_llk +EGGllk +H_ADFGllk +H_PWSSCllk +MDMllk +age0_devs_penllk +mort_devs_penllk +age0_covar_prior +mort_covar_prior +Z_prior +hydADFG_add_prior +hydPWSSC_add_prior +m_add_prior +juv_llk +sero_llk;

  SSB_final_year = SSB(nyr_tobefit); // is projected year's Pre-fishery Run Biomass in metric tons
  
  // "ifMCEvalPhase" goes inside the PROCEDURE_SECTION,
  if(mceval_phase()){
    project_biomass(); // Project current year biomass using last year's data
    write_chain_results();
  }

// |---------------------------------------------------------------------------|


FUNCTION void runSimulationModel(const int& rseed)
  //PSEUDOCODE
  //  1) Calculate mortality, maturity, and selectivity (independent of biomass)
  //  2) If selected, calculate catches from user-input exploitation history.
  //    i) Selectivities of fisheries are approximated on biomass calculating state variables with actual historical data & estimated Age 0 devs (not random) 
  //  3) If selected, randomly generate Age 0 devs
  //  4) If selected, calculate averages of data not fit in the model (weight-at-age, fecundity-at-age, % female)
  //  5) Calculate simulated state variables (biomass)
  //  6) Calculate survey data (unobserved values)
  //  7) Simulate observed survey data based on historical SE's or ESS
  //  8) Output data files with new data AND true biomass values

  Mean_Age0 = exp(log_MeanAge0+Mean_Age0offset*R_change(1,nyr_tobefit));
  mdm_c = exp(logmdm_c);
  init_pop = exp(loginit_pop);

  calc_naturalmortality();
  if(DD_Mat==0){
	  calc_maturity();
  }
  calc_selectivity();
  
  random_number_generator rng(rseed);

  // Following allows for conditioning on user-specified target harvest rate instead of directly from catches
  // This can be specified from command line
  if(sim_catches>0){

    calc_statevariables();
    // cout<<"FOR ERROR CHECKING"<<endl<<endl;
    // We want to include 0 in case simulations are in absence of fishing

    // Approach is somewhat non-parametric simulation of catch-at-age 
    // Use historical catch-at-age and estimated Numbers at age from projections with
    // historical catches to calculate proxy selectivities for each year

    // One row for each fishery
    dmatrix prop_harvest(1,4,1,nyr_tobefit); 
    // dvar_matrix age_comp(1,nyr,1,nage);
    // dvar_vector total_historical(1,nyr_tobefit);

    prop_harvest(1)=exploitation_history/4;  // For now assume that harvest is evenly divied up to the 4 fisheries
    prop_harvest(2)=exploitation_history/4; 
    prop_harvest(3)=exploitation_history/4; 
    prop_harvest(4)=exploitation_history/4; 

    // Gillnet fishery vulnerability
    // total_historical =rowsum(gc.sub(1,nyr_tobefit))+0.000000001; // To prevent errors from dividing by 0 for no catch years
    for(int i=1; i<=nyr_tobefit; i++){
      // age_comp(i)(1,nage) = gc(i)(1,nage)/total_historical(i);

      // Derived from equation used to calculate purse seine age comps for which only total catch and estimated selectivity is provided
      // gc_V(i)(1,nage) = value(elem_div(age_comp(i)(1,nage)*total_historical(i),N_y_a(i)(1,nage)+0.000000001)); 

      // Multiply N_y_a again in numerator and denominator for those estimated 0 years in upper age classes at the beginning of the time series
      gc_V(i)(1,nage) = value(elem_div(elem_prod(gc(i)(1,nage),N_y_a(i)(1,nage)),elem_prod(N_y_a(i)(1,nage),N_y_a(i)(1,nage))+0.000000001)); 

      gc_V(i)(1,nage) = gc_V(i)(1,nage)/max(gc_V(i)(1,nage));

      // Assumptions about this approach: resulting curve approximates annual selectivity & maximum value represents fully selected age in the fishery
    }

    // Resample from only those years fished to fill gc_V in all years - allow user to define time frame from which to resample
    // sample(vector of samples,# of samples in output,1 with replacement or 0 without replacement,rng)
    ivector resampled_years(1,nyr_tobefit);
    resampled_years=resample_period(sample(resample_period,nyr_tobefit,1,rng));
    
    for(int i=1; i<=nyr_tobefit; i++){
      gc_V(i)(1,nage) = prop_harvest(1,i)*gc_V(resampled_years(i))(1,nage);
    }


    // Impound fishery vulnerability
    for(int i=1; i<=nyr_tobefit; i++){
      pc_V(i)(1,nage) = value(elem_div(elem_prod(pc(i)(1,nage),N_y_a(i)(1,nage)),elem_prod(N_y_a(i)(1,nage),N_y_a(i)(1,nage))+0.000000001)); 
      pc_V(i)(1,nage) = pc_V(i)(1,nage)/max(pc_V(i)(1,nage));
    }
    for(int i=1; i<=nyr_tobefit; i++){
      pc_V(i)(1,nage) = prop_harvest(2,i)*pc_V(resampled_years(i))(1,nage);
    }

    // Food & Bait fishery vulnerability
    dvar_vector temp_N_y_a(1,nage);
    for(int i=1; i<=nyr_tobefit; i++){
      // Need to account for other catches since this fishery is in 2nd half of model year
      temp_N_y_a = elem_prod(N_y_a(i)(1,nage)-(SeAC(i)(1,nage)*N_se(i)+gc(i)(1,nage)+pk*pc(i)(1,nage)),Sur_summer(i)(1,nage));
      // Multiply numerator and denominator by N_y_a again to account for beginning years where N_y_a is 0 (then temp_N_y_a is negative)
      fbc_V(i)(1,nage) = value(elem_div(elem_prod(fbc(i)(1,nage),N_y_a(i)(1,nage)),elem_prod(temp_N_y_a,N_y_a(i)(1,nage))+0.000000001)); 
      fbc_V(i)(1,nage) = fbc_V(i)(1,nage)/max(fbc_V(i)(1,nage));
    }

    for(int i=1; i<=nyr_tobefit; i++){
      fbc_V(i)(1,nage) = prop_harvest(3,i)*fbc_V(resampled_years(i))(1,nage);
    }

    // Seine fishery - since selectivity is estimated (i.e. a parameter) in the model, just need annual seine yield
    sc_F = prop_harvest(4)(1,nyr_tobefit);
    
    turn_on_effort = 1;

    // Notes for future development
    // Use dirichlet to generate catches
    // rdirichlet(shape,rng); returns a vector that is rgamma(shape[i],rng)[i]/sum of all rgamma 
    // shape vector is age distribution, either as proportions or actual numbers
    // Pseudocode for simulated catches
    // 3) Two other options to include:
    //    b) Annual age comps divided by N_y_a, scale each age by maximum, calculate average selectivity for period of years and fill Nyr x Nage matrix (can have option to allow dirichlet draws of age comps)
    //    d) Assume functional selectivity - logistic or dome-shaped selectivity, and carry out calculations like the purse-seine age comp calculations (w/ dirichlet)
  }


  // if age0_dev_option set to 1, then random process error added
  // Otherwise, devs from PIN file are used
  if(age0_dev_option==1){
    dvector simulated_age0devs(1,nyr_tobefit-3);
    simulated_age0devs.fill_randn(rng);
    annual_age0devs(1)(1,nyr_tobefit-3) = dvar_vector(simulated_age0devs * sigma_age0devs-0.5*square(sigma_age0devs));
  }


  // Several data sets are not fit by the model - weight & fecundity
  // This option allows to use the average over all years of these data for simulations
  if(data_avg_option>0){
    dvar_vector temp_fecun(1,nage);
    int temp_yrcount=0;

    temp_fecun=0;
    for(int i=1; i<=nyr_tobefit; i++){
      if(sum(fecun(i)(1,nage))>0){
        temp_fecun+=fecun(i)(1,nage);
        temp_yrcount+=1;
      }
    }

    dvar_vector temp_waa(1,nage);
    double temp_f_sp=sum(f_sp(1,nyr_tobefit))/nyr_tobefit;
    
    temp_waa=colsum(w_a_a.sub(1,nyr_tobefit))/nyr_tobefit;
    for(int i=1; i<=nyr_tobefit; i++){
      w_a_a(i)(1,nage)=value(temp_waa);
      fecun.rowfill(i,value(temp_fecun/temp_yrcount));
      f_sp(i) = temp_f_sp;
    }
  }

  calc_statevariables();
  calc_surveyvalues();

  // Observation model - rewrite data inputs, change to simulated data
  // If no rseed is provided, default is to store unobserved (true) values 
  // 12/2019 - need to incorporate simulation of juvenile survey index
  if(rseed>0){
    dvector egg_obs_error(1,nyr_tobefit);
    egg_obs_error.fill_randn(rng);

    dvector milt_obs_error(1,nyr_tobefit);
    milt_obs_error.fill_randn(rng);

    dvector ADFG_acoustic_obs_error(1,nyr_tobefit);
    ADFG_acoustic_obs_error.fill_randn(rng);

    dvector PWSSC_acoustic_obs_error(1,nyr_tobefit);
    PWSSC_acoustic_obs_error.fill_randn(rng);

    for(int i=1; i<=nyr_tobefit; i++){
      if(egg(i)>0){egg(i) = value(EGG(i)*exp(egg_obs_error(i)*Eg_SD(i)-0.5*square(Eg_SD(i))));}

      if(mdm(i)>0){mdm(i) = value(MDM(i)*exp(milt_obs_error(i)*m_add-0.5*square(m_add)));}

      if(hydADFG(i)>0){hydADFG(i) = value(HYD_ADFG(i)*exp(ADFG_acoustic_obs_error(i)*hydADFG_add-0.5*square(hydADFG_add)));}

      if(hydPWSSC(i)>0){hydPWSSC(i) = value(HYD_PWSSC(i)*exp(PWSSC_acoustic_obs_error(i)*PWSSC_SD(i)-0.5*square(PWSSC_SD(i))));}

      // Rewrite age comps
      if(seine(i,4)>=0){seine(i)(4,nage) = rmultinom(rseed+i,ESS_Se(i),value(SeAC(i)(4,nage)))/ESS_Se(i);}
      if(spac(i,4)>=0){spac(i)(4,nage) = rmultinom(rseed+i,ESS_Sp(i),value(SpAC(i)(4,nage)))/ESS_Sp(i);}
    }
  }else{
    for(int i=1; i<=nyr_tobefit; i++){
      if(egg(i)>0){egg(i) = value(EGG(i));}
      if(mdm(i)>0){mdm(i) = value(MDM(i));}
      if(hydADFG(i)>0){hydADFG(i) = value(HYD_ADFG(i));}
      if(hydPWSSC(i)>0){hydPWSSC(i) = value(HYD_PWSSC(i));}

      if(seine(i,4)>=0){seine(i)(4,nage) = value(SeAC(i)(4,nage));}
      if(spac(i,4)>=0){spac(i)(4,nage) = value(SpAC(i)(4,nage));}
    }
  }

  ofstream sim_data("PWS_ASA_sim.dat",ios::trunc);
  ofstream sim_state_vars("simulated_biomass.dat",ios::trunc);
  
  sim_state_vars << "# Pre-fishery Spawning Biomass (metric tons)" << endl;
  sim_state_vars << SSB << endl << endl;
  sim_state_vars << "# Post-fishery Spawning Biomass (metric tons)" << endl;
  sim_state_vars << SB << endl << endl;
  sim_state_vars << "# Recruitment (millions of Age 3 herring)" << endl;
  for(int i=1; i<=nyr_tobefit; i++){
    sim_state_vars << N_y_a(i,4) << " ";
  }
  sim_state_vars << endl;

  // 12/2019 - Will need to edit to exactly match PWS_ASA.dat file (there are other things I added)
  sim_data << "# Simulated data from PWS BASA model" << endl << endl;
  sim_data << "# Number of years (nyr)" << endl;
  sim_data << nyr_tobefit << endl << endl;
  sim_data << "# Number of years to be fit in the model (nyr_tobefit)" << endl;
  sim_data << nyr_tobefit << endl << endl;
  sim_data << "#Number of age classes (nage)" << endl;
  sim_data << nage << endl << endl;
  sim_data << "#Weight-at-age (w_a_a)" << endl;
  sim_data << w_a_a.sub(1,nyr_tobefit) << endl << endl;
  sim_data << "#Fecundity-at-age (fecun)" << endl;
  sim_data << fecun.sub(1,nyr_tobefit) << endl << endl;
  sim_data << "# Pound Catch (pc)" << endl;
  sim_data << pc.sub(1,nyr_tobefit) << endl << endl;
  sim_data << "#Proportion of pound fish killed (pk)" << endl;
  sim_data << pk << endl << endl;
  sim_data << "#Food/bait catch (fbc)" << endl;
  sim_data << fbc.sub(1,nyr_tobefit) << endl << endl;
  sim_data << "#Gillnet catch (gc)" << endl;
  sim_data << gc.sub(1,nyr_tobefit) << endl << endl;
  sim_data << "#Total Seine yield (sc)" << endl;
  sim_data << sc(1,nyr_tobefit) << endl << endl;
  sim_data << "# % Female (fem)" << endl;
  sim_data << f_sp(1,nyr_tobefit) << endl << endl;
  sim_data << "# Mile-days of mile (mdm)" << endl;
  sim_data << mdm(1,nyr_tobefit) << endl << endl;
  sim_data << "# Egg deposition (egg)" << endl;
  sim_data << egg(1,nyr_tobefit) << endl << endl;
  sim_data << "# Egg Deposition standard errors - only 10 years - assuming C.I.'s used Normal" << endl;
  sim_data << cv_egg(1,nyr_tobefit) << endl << endl;
  sim_data << "# Index of first year of hydroacoustic survey from ADFG" << endl;
  sim_data << hydADFG_start << endl << endl;
  sim_data << "# Hydroacoustic survey (hyd) from ADFG" << endl;
  sim_data << hydADFG(1,nyr_tobefit) << endl << endl;
  sim_data << "# Index of first year of hydroAcoustic survey for PWSSC" << endl;
  sim_data << hydPWSSC_start << endl << endl;
  sim_data << "# Hydroacoustic survey (hyd) from PWSSC " << endl;
  sim_data << hydPWSSC(1,nyr_tobefit) << endl << endl;
  sim_data << "# PWSSC Hydroacoustic Biomass s.e." << endl;
  sim_data << cv_hydPWSSC(1,nyr_tobefit) << endl << endl;
  sim_data << "#Seine age composition (seine)" << endl;
  sim_data << seine.sub(1,nyr_tobefit) << endl << endl;
  sim_data << "#Spawning age composition (spac)" << endl;
  sim_data << spac.sub(1,nyr_tobefit) << endl << endl;



FUNCTION void calc_naturalmortality()
  summer_effect.initialize();
  winter_effect.initialize();
  Sur_summer.initialize();
  Sur_winter.initialize();
  annual_mortdevs_byage.initialize();
   
  //Half-Year Survival (Matches Excel version - uses desease data and only estimates plus group mortality)
  //S_3_8=exp(-0.5*Z_0_8); //Z_0_8 is a read-in param 
  
  // below should be added to sdreport? Nah, but should look at aspects the survival matrix as sdreport candidates.
  S_9=exp(-0.5*Z_9); // Z_9 is estimated
 
// This sets up the survival matrix
// I only go up to i<nyr_tobefit because the model's calendar is different from the normal one
// Every new model year starts at time of spawn (spring), so first half-year survival is summer,
// and second half-year survival is winter. This lag of a year is accounted for in the next bit.
// if we are predicting herring in nyr_tobefit, we consider any summer mortality effect from nyr_tobefit-1 and 
// winter mortality effect occurs in this nyr_tobefit (the same year). But, for indexing survival with i, 
// survival for year i is the summer mortality in calendar year i and winter mortality in calendar
// year i+1
  
  for (int i=1;i<=nyr_tobefit;i++){
    // Z_annual(i)=(M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(k)*mor_covariates(i,k)*covariate_effect_byage(j,k);
    for(int j=1;j<=nage;j++){
      for(int k=1;k<=mor_cov_counter;k++){
        if(beta_mortality_ind(k)==0){
        }else if(mor_covariates(i,beta_mortality_ind(k))==-9){
          summer_effect(i,j) += 0;
          winter_effect(i,j) += 0;
        }else if(mor_season(beta_mortality_ind(k))==1){
          summer_effect(i,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
        }else if(mor_season(beta_mortality_ind(k))==2){
          winter_effect(i,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
        }else if(mor_season(beta_mortality_ind(k))==3){
        // Additional conditional because each model year starts with summer, ends with winter even though 
        // data are input for the calendar year
          if(i==nyr_tobefit){
            summer_effect(i,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
          }else{
            summer_effect(i,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
            winter_effect(i+1,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
          } 
        }
      }
    }
  }

// Needed because calendar year for winter effect is in the second year-season of the model. In other words, user can input 
// covariate effect info for the year being forecasted.

  for(int j=1;j<=nage;j++){
    for(int k=1;k<=mor_cov_counter;k++){
      if(beta_mortality_ind(k)==0){
      }else if(nyr_tobefit_winter_covariate(k)==-9){
        forecast_winter_effect(j) += 0;
      }else if(mor_season(beta_mortality_ind(k))==2){
        forecast_winter_effect(j) += (M_change(nyr_tobefit)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*nyr_tobefit_winter_covariate(beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
      }else if(mor_season(beta_mortality_ind(k))==3){
        forecast_winter_effect(j) += (M_change(nyr_tobefit)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*nyr_tobefit_winter_covariate(beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
      }
    }
  }

  // Added 10/16/2018: Find which ages are not affected by beta - FIX to Z_0_8 (assumed to not be estimated while Z_0_8_offset is estimated to center with beta estimate)
  dvar_vector Fix_Z_summer(1,nage);
  dvar_vector Fix_Z_winter(1,nage);
  Fix_Z_summer.initialize();
  Fix_Z_winter.initialize();
  for(int j=1;j<=nage;j++){
    for(int k=1;k<=n_mor_covs;k++){
      if(covariate_effect_byage(j,k)*mor_turn_on(k)==1){
      	if(Fix_Z_summer(j)!=1 & (mor_season(k)==1 | mor_season(k)==3)){
      		Fix_Z_summer(j)=1;
      	}
      	if(Fix_Z_winter(j)!=1 & (mor_season(k)==2 | mor_season(k)==3)){
      		Fix_Z_winter(j)=1;
      	}
  	  }
  	}
  }
  
  
  if(M_cov_model==1){
    // Previous form where covariates are incoporated as fixed variables - changed 07/05/2019
    for (int i=1;i<=nyr_tobefit;i++){
      for (int j=1;j<=(nage-1);j++){
        if(i==13){
          if((j>=4) && (j<=5)){
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)+VHSV_age3_4_mort_93));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)));

          }else if(j>=6){
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)+ICH_age5_8_mort_93));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)));
          }else{
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)));
          }
        }else if(i==14){
          if((j>=4) && (j<=5)){
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)+VHSV_age3_4_mort_93));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)+VHSV_age3_4_mort_93));
          }else if(j>=6){
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)+ICH_age5_8_mort_93));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)+ICH_age5_8_mort_93));
          }else{
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)));
          }
        }else if(i==15){
          if((j>=4) && (j<=5)){
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)+VHSV_age3_4_mort_93));
          }else if(j>=6){
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)+ICH_age5_8_mort_93));
          }else{
            Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)));
            Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)));
          }
        }else{
          Sur_summer(i,j)=exp(-((0.5*M_change(i)*Fix_Z_summer(j)*Z_0_8offset+0.5*Z_0_8)+summer_effect(i,j)));
          Sur_winter(i,j)=exp(-((0.5*M_change(i)*Fix_Z_winter(j)*Z_0_8offset+0.5*Z_0_8)+winter_effect(i,j)));
        }
        
        // Calculate survival for ages 3-8 for years 1981 and above  
        dvariable pen_Sur_1=0.0;
        dvariable high_survival_penalty_1=1-Sur_summer(i,j);
        Sur_summer(i,j)=1-posfun(high_survival_penalty_1, 0.01, pen_Sur_1);
          if(Sur_summer(i,j)>=1){
            penCount+=1;
          }
          f_llk+=1000*pen_Sur_1;

        dvariable pen_Sur_2=0.0;
        dvariable high_survival_penalty_2=1-Sur_winter(i,j);
        Sur_winter(i,j)=1-posfun(high_survival_penalty_2, 0.01, pen_Sur_2);
          if(Sur_winter(i,j)>=1){
            penCount+=1;
          }
          f_llk+=1000*pen_Sur_2;
      }

      if(i==1){
        Sur_summer(i,nage)=exp(-((0.5*M_change(i)*Fix_Z_summer(nage)*Z_9offset+0.5*Z_9)+summer_effect(i,nage)));
        Sur_winter(i,nage)=exp(-((0.5*M_change(i)*Fix_Z_winter(nage)*Z_9offset+0.5*Z_9)+winter_effect(i,nage)));

        dvariable pen_Sur_3=0.0;
        dvariable high_survival_penalty_3=1-Sur_summer(i,nage);
        Sur_summer(i,nage)=1-posfun(high_survival_penalty_3, 0.01, pen_Sur_3);
          if(Sur_summer(i,nage)>=1){
            penCount+=1;
          }
          f_llk+=1000*pen_Sur_3;

        dvariable pen_Sur_4=0.0;
        dvariable high_survival_penalty_4=1-Sur_winter(i,nage);
        Sur_winter(i,nage)=1-posfun(high_survival_penalty_4, 0.01, pen_Sur_4);
         if(Sur_winter(i,nage)>=1){
           penCount+=1;
         }
        f_llk+=1000*pen_Sur_4;
      }else {
        Sur_summer(i,nage)=Sur_summer(i,nage-1)*Sur_summer(i-1,nage)/Sur_summer(i-1,nage-1); // Plus age group
        Sur_winter(i,nage)=Sur_winter(i,nage-1)*Sur_winter(i-1,nage)/Sur_winter(i-1,nage-1); // Plus age group
      }
    }

  }else if(M_cov_model==2){
    for (int i=1;i<=nyr_tobefit;i++){
      for (int j=1;j<=(nage-1);j++){
        // I must do this because I only want to estimate deviates on mortality for ages and years to which I am fitting data
        for (int k=1;k<=mor_cov_counter;k++){
          if(beta_mortality_ind(k)==0){
          }else if(mor_covariates(i,beta_mortality_ind(k))!=-9){
            annual_mortdevs_byage(i,j)+=(M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*covariate_effect_byage(j,beta_mortality_ind(k))*annual_mortdevs(k,i);
          }
        }
        Sur_summer(i,j)=exp(-(0.5*(M_change(i)*Fix_Z_summer(j)*Z_0_8offset+Z_0_8)+annual_mortdevs_byage(i,j)));
        Sur_winter(i,j)=exp(-(0.5*(M_change(i)*Fix_Z_winter(j)*Z_0_8offset+Z_0_8)+annual_mortdevs_byage(i,j)));

        // Calculate survival for ages 3-8 for years 1981 and above  
        dvariable pen_Sur_1=0.0;
        dvariable high_survival_penalty_1=1-Sur_summer(i,j);
        Sur_summer(i,j)=1-posfun(high_survival_penalty_1, 0.01, pen_Sur_1);
          if(Sur_summer(i,j)>=1){
            penCount+=1;
          }
          f_llk+=1000*pen_Sur_1;

        dvariable pen_Sur_2=0.0;
        dvariable high_survival_penalty_2=1-Sur_winter(i,j);
        Sur_winter(i,j)=1-posfun(high_survival_penalty_2, 0.01, pen_Sur_2);
          if(Sur_winter(i,j)>=1){
            penCount+=1;
          }
          f_llk+=1000*pen_Sur_2;
      }

      if(i==1){

        Sur_summer(i,nage)=exp(-(0.5*(M_change(i)*Fix_Z_summer(nage)*Z_9offset+Z_9)+annual_mortdevs_byage(i,nage)));
        Sur_winter(i,nage)=exp(-(0.5*(M_change(i)*Fix_Z_winter(nage)*Z_9offset+Z_9)+annual_mortdevs_byage(i,nage)));

        dvariable pen_Sur_3=0.0;
        dvariable high_survival_penalty_3=1-Sur_summer(i,nage);
        Sur_summer(i,nage)=1-posfun(high_survival_penalty_3, 0.01, pen_Sur_3);
          if(Sur_summer(i,nage)>=1){
            penCount+=1;
          }
          f_llk+=1000*pen_Sur_3;

        dvariable pen_Sur_4=0.0;
        dvariable high_survival_penalty_4=1-Sur_winter(i,nage);
        Sur_winter(i,nage)=1-posfun(high_survival_penalty_4, 0.01, pen_Sur_4);
         if(Sur_winter(i,nage)>=1){
           penCount+=1;
         }
        f_llk+=1000*pen_Sur_4;
      }else {
        Sur_summer(i,nage)=Sur_summer(i,nage-1)*Sur_summer(i-1,nage)/Sur_summer(i-1,nage-1); // Plus age group
        Sur_winter(i,nage)=Sur_winter(i,nage-1)*Sur_winter(i-1,nage)/Sur_winter(i-1,nage-1); // Plus age group
      }
    }     
  }
  


FUNCTION void calc_maturity()
  Mat.initialize();
  Mat_unobs.initialize();
  
  if(mat_mod_type==1){
    // Maturity values before 1997
    // Change to 14 if pre 1994, 17 pre 1997, nyr_tobefit if single Early_biomass
    int yr_mat_change=nyr_tobefit;
    // int yr_mat_change=17;

    for (int i=1;i<=nyr_tobefit;i++){
      Mat(i)(1,3)=0; 
      if(i<=yr_mat_change){
        // Mat(i,4)=matur_age3_per1*matur_age4_per1;
        // Mat(i,5)=matur_age4_per1;

        // Logistic maturity function
        // matur_age3_per1 = age @ 50% mature
        // matur_age4_per1 = age @ 95% mature
        Mat(i,4)=1/(1+exp(-log(19)*(3-matur_age3_per1)/(matur_age4_per1-matur_age3_per1)));
        Mat(i,5)=1/(1+exp(-log(19)*(4-matur_age3_per1)/(matur_age4_per1-matur_age3_per1)));
      }else {
        // Mat(i,4)=matur_age3_per2;
        // Mat(i,5)=matur_age4_per2;

        Mat(i,4)=1/(1+exp(-log(19)*(3-matur_age3_per2)/(matur_age4_per2-matur_age3_per2)));
        Mat(i,5)=1/(1+exp(-log(19)*(4-matur_age3_per2)/(matur_age4_per2-matur_age3_per2)));
      }
      Mat(i,6)=1;
      Mat(i)(7,nage)=1;
    }
  }else{
    for (int i=1;i<=nyr_tobefit;i++){
      Mat(i)(1,3)=0; 
      Mat(i,4)=matur_age3_per1*matur_age4_per1;
      Mat(i,5)=matur_age4_per1;
      Mat(i)(6,nage)=1;

      Mat_unobs(i)(1,3)=0; 
      Mat_unobs(i,4)=matur_age3_per2*matur_age4_per2;
      Mat_unobs(i,5)=matur_age4_per2;
      Mat_unobs(i)(6,nage)=1;
    }
  }
  


FUNCTION void calc_selectivity()
  
  //Gear Selectivity
  Vul.initialize();
  Vul(1,3)=0; 
  for (int j=4;j<=nage;j++) {
    Vul(j)=1/(1+exp(-1.0*beta_v*(j-1-alpha_v)));
  }

  //Survey Selectivity
  Vul_survey.initialize();
  Vul_survey(1,3)=0; 
  for (int j=4;j<=nage;j++) {
    Vul_survey(j)=1/(1+exp(-1.0*survey_vul_beta*(j-1-survey_vul_alpha)));
  }

FUNCTION void calc_statevariables()
  N_y_a.initialize();
  CC_top.initialize();
  CC_bot.initialize();
  SeAC.initialize();
  bio_a_a.initialize();
  SeBiomass.initialize();
  N_se.initialize();
  N_sp.initialize();
  Early_Sp_biomass.initialize();
  SSB.initialize();
  SB_star.initialize();
  SPB.initialize();
  SB.initialize();
  SpAC_top.initialize();
  SpAC_bot.initialize();
  SpAC.initialize();
  Early_bio_a_a.initialize();
  Early_biomass.initialize();
  age0_effect.initialize();
  init_age_0.initialize();
  agggregate_annual_age0devs.initialize();
  suscep.initialize();
  immune.initialize();

  not_below_this = 0.01;

  // Age-specific vulnerability to sympatric transmission
  Vul_sympat.initialize();
  for (int a=1;a<=nage;a++){
    // Vul_sympat(a) = 1/(1+exp(-log(19)*((a-1)-infec_vul_a50)/(infec_vul_a95-infec_vul_a50)));
    Vul_sympat(a) = 1.0;
  }

  // All new born fish are 100% susceptible - COME BACK
  // suscep.colfill(1,1.0);
  suscep = 1.0;
  Sur_vhsv = 1.0;

  // Time-varying OR constant recovery probability
  dvar_vector pseudo_rec(1,nyr_tobefit);
  if(recov_prob_type==1){
    pseudo_rec = recov_prob(1);  // When assuming constant, only uses first parameter in vector
  }else{
    pseudo_rec = recov_prob;
  }

  // INITIAL CONDITIONS
  // init_age_0 is the number of recruits in year i
    // Sur_age0_2 is the survival rate experienced by init_age_0(i) from age 0.
    // In other words, the mortality the larvae in year i had undergone. Referred to as recruit because previous model started at age 3
    for(int k=1;k<=rec_cov_counter;k++){
      if(beta_recruit_ind(k)==0){
      }else if(age0_covariates(1,beta_recruit_ind(k))==-9){
        age0_effect(1) += 0;
      }else{
        age0_effect(1) += age0_turn_on(beta_recruit_ind(k))*(R_change(1)*beta_age0_offset(k)+beta_age0(k))*age0_covariates(1,beta_recruit_ind(k));
        
        // This is only used if R_cov_model==2
        agggregate_annual_age0devs(1) += (R_change(1)*beta_age0_offset(k)+beta_age0(k))*annual_age0devs(k,1);
      }
    }

    if(R_cov_model==1){
      // Form below for effects incorporated as fixed variables
      init_age_0(1) = Mean_Age0(1)*exp((age0_effect(1)+annual_age0devs(1,1)-0.5*square(sigma_age0devs))); 
    }else if(R_cov_model==2){
      init_age_0(1) = Mean_Age0(1)*exp(agggregate_annual_age0devs(1)-0.5*square(sigma_age0devs)); 
    }

  //Fills in row 1 of pre-fishery abundance
    N_y_a(1,1)=init_age_0(1);
    --N_y_a(1)(2,6)=init_pop;   // Fill in first year N_y_a subvector with estimated initial ages 1-5 
    N_y_a(1)(7,nage)=0;      

  if(DD_Mat==1){
  	Mat.initialize();
  	Mat(1)(1,3)=0;
  	Mat(1,4)=1/(1+exp(matur_age3_per1+exp(matur_age4_per1)*sum(N_y_a(1)(1,nage))));
  	Mat(1,5)=Mat(1,4)+(1-Mat(1,4))/(1+exp(matur_age3_per2));
  	Mat(1)(6,nage)=1;
  }
  // Calculate pre-fishery spawning biomass - nest elem_prod because function only accepts 2 args
  if(mat_mod_type==1){
    Early_Sp_biomass(1)(1,nage)=elem_prod(elem_prod(Mat(1)(1,nage),N_y_a(1)(1,nage)),w_a_a(1)(1,nage));
  }else{
    for(int j=1;j<=nage;j++){
      Early_Sp_biomass(1,j)=(Vul_survey(j)*Mat(1,j)+(1-Vul_survey(j))*Mat_unobs(1,j))*N_y_a(1,j)*w_a_a(1,j);
    }
  }

  
  SSB(1)=sum(Early_Sp_biomass(1)(1,nage));

  SB_star[1]=SSB[1]*2204.62/2000; // pre-fishery run biomass in TONS
  
  //Fills in row 1 of Catch Age-Composition from Purse-Seine
  CC_top(1)(1,nage)=elem_prod(Vul,N_y_a(1)(1,nage));
  CC_bot(1)=sum(CC_top(1)(1,nage)); 

  SeAC(1)(1,nage)=CC_top(1)(1,nage)/CC_bot(1);
  bio_a_a(1)(1,nage)=elem_prod(SeAC(1)(1,nage),w_a_a(1)(1,nage));
  SeBiomass=sum(bio_a_a(1));

  if(mat_mod_type==1){
    SpAC_top(1)(1,nage)=elem_prod(Mat(1)(1,nage),N_y_a(1)(1,nage)); //numerator of Spawning Age-Comp(SpAC)
  }else{
    SpAC_top(1)(1,nage)=elem_prod(Vul_survey(1,nage),N_y_a(1)(1,nage)); //numerator of Spawning Age-Comp(SpAC)
  }
  SpAC_bot(1)=sum(SpAC_top(1)(1,nage));
  SpAC(1)(1,nage)=SpAC_top(1)(1,nage)/SpAC_bot(1);  //fills in Spawning Age-Comp

  if(turn_on_effort){
    sc = 0;// Set all catches to 0 because we calculate them as we go along
    // For simulation, calculate new seine yield based on specified exploitation/harvest rate sc_F from simulation section
    sc(1)=value(sc_F(1)*sum(elem_prod(Vul(1,nage),elem_prod(N_y_a(1)(1,nage),w_a_a(1)(1,nage)))));
  }
  //Row 1 of total Seine catch in millions 
  N_se(1)=sc(1)/SeBiomass(1);

  int count = 7;
  double gc_sum;
  double pc_sum;
  double fbc_sum;
    
  // Now generate naturally spawning pop as correct structure for first four years
  N_sp(1)(1,nage)=0;
  SPB(1)(1,nage)=0;

  SB(1)=0;

  if(turn_on_effort){
    gc = 0;
    pc = 0;
    fbc = 0; // Set all catches to 0 because we calculate them as we go along

    gc(1)(1,nage)=value(elem_prod(gc_V(1)(1,nage),N_y_a(1)(1,nage)));
    pc(1)(1,nage)=value(elem_prod(pc_V(1)(1,nage),N_y_a(1)(1,nage)));
    fbc(1)(1,nage)=value(elem_prod(fbc_V(1)(1,nage),elem_prod(N_y_a(1)(1,nage)-(SeAC(1)(1,nage)*N_se(1)+gc(1)(1,nage)+pk*pc(1)(1,nage)),Sur_summer(1)(1,nage))));
  }

  N_sp(1)(4,6)=elem_prod(Mat(1)(4,6),N_y_a(1)(4,6)-(SeAC(1)(4,6)*N_se(1)+gc(1)(4,6)+pk*pc(1)(4,6)));

  for(int j=4;j<=6;j++){
    dvariable pen4=0.0;
    N_sp(1,j)=posfun(N_sp(1,j), .01, pen4);
    if(N_sp(1,j) <= not_below_this){
      penCount+=1;
    }
    f_llk+=1000*pen4;
  } 
  
  //SPB(1)(4,6)=elem_prod(elem_prod(N_sp(1)(4,6),w_a_a(1)(4,6)),1/Mat(1)(4,6));
  SPB(1)(4,6)=elem_prod(N_sp(1)(4,6),w_a_a(1)(4,6));
  
  SB(1)=sum(SPB(1)(4,6)); // SB=rowsum(SPB); //Total naturally spawning biomass
  int m=7;

  // FILL IN THE REMAINING YEARS

  for(int i=2;i<=nyr_tobefit;i++){
    // init_age_0 is the number of recruits in year i

    if(i<=(nyr_tobefit-3)){
      temp_age0 = column(annual_age0devs,i);
    }else{
      temp_age0 = 0;
    }

    for(int k=1;k<=rec_cov_counter_age0devs;k++){
      if(beta_recruit_ind(k)==0){
      }else if(age0_covariates(i,beta_recruit_ind(k))==-9){
        age0_effect(i) += 0;
      }else{
        age0_effect(i) += age0_turn_on(beta_recruit_ind(k))*(R_change(i)*beta_age0_offset(k)+beta_age0(k))*age0_covariates(i,beta_recruit_ind(k));

        agggregate_annual_age0devs(i) += (R_change(i)*beta_age0_offset(k)+beta_age0(k))*temp_age0(k);
      }
    }

    if(R_cov_model==1){
      init_age_0(i) = Mean_Age0(i)*exp((age0_effect(i)+temp_age0(1)-0.5*square((sigma_age0devs))));
    }else if(R_cov_model==2){
      init_age_0(i) = Mean_Age0(i)*exp(agggregate_annual_age0devs(i)-0.5*square(sigma_age0devs)); 
    }

    if(i<=5){                     //Fill in years 2:5 as plus group advances from 6+ to 9+
      N_y_a(i,1)=init_age_0(i);
      for(int j=2;j<=count-1;j++){
        if(turn_on_effort){
          gc(i-1,j-1)=value(gc_V(i-1,j-1)*N_y_a(i-1,j-1));
          pc(i-1,j-1)=value(pc_V(i-1,j-1)*N_y_a(i-1,j-1));
          fbc(i-1,j-1)=value(fbc_V(i-1,j-1)*(N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc(i-1,j-1)+pk*pc(i-1,j-1)))*Sur_summer(i-1,j-1));
        }
        N_y_a(i,j)=((N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc(i-1,j-1)+pk*pc(i-1,j-1)))*Sur_summer(i-1,j-1)-fbc(i-1,j-1))*Sur_winter(i,j-1);

      }
      //The last age-class each year (incrementally by year beginning with age 6 in 1981) is a plus group
      gc_sum = gc(i-1,count-1);
      pc_sum = pc(i-1,count-1);
      fbc_sum = fbc(i-1,count-1);
      for(int k=count;k<=nage;k++){
        if(turn_on_effort){
          gc(i-1,k)=value(gc_V(i-1,k)*N_y_a(i-1,k));
          pc(i-1,k)=value(pc_V(i-1,k)*N_y_a(i-1,k));
          fbc(i-1,k)=value(fbc_V(i-1,k)*(N_y_a(i-1,k)-(SeAC(i-1,k)*N_se(i-1)+gc(i-1,k)+pk*pc(i-1,k)))*Sur_summer(i-1,k));
        }
        gc_sum += gc(i-1,k);
        pc_sum += pc(i-1,k);  
        fbc_sum += fbc(i-1,k); 
      }
      int j=count;

      // Ignore disease survival here because no disease data
      N_y_a(i,j)=(((N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc_sum+pk*pc_sum))*Sur_summer(i-1,j-1))-fbc_sum)*Sur_winter(i,j-1); 

      for(int j=1;j<=count;j++){
        dvariable pen1=0.0;
        N_y_a(i,j)=posfun(N_y_a(i,j), 0.01, pen1);
        if(N_y_a(i,j)<=not_below_this){
          penCount+=1;
        }       
       f_llk+=1000*pen1;
       CC_top(i,j)=N_y_a(i,j)*Vul(j);
      }

      CC_bot=rowsum(CC_top); 
      for(int j=1;j<=count;j++){
         SeAC(i,j)=CC_top(i,j)/CC_bot(i);
         bio_a_a(i,j)=SeAC(i,j)*w_a_a(i,j); 
      }
      SeBiomass=rowsum(bio_a_a);
      
      if(turn_on_effort){
        // For simulation, calculate new seine yield based on specified exploitation/harvest rate sc_F from simulation section
        sc(i)=value(sc_F(i)*sum(elem_prod(Vul(1,nage),elem_prod(N_y_a(i)(1,nage),w_a_a(i)(1,nage)))));
      }
                
      N_se(i)=sc(i)/SeBiomass(i);
      count+=1;

    } else if((i>5) && (i<=13)){
      //Fills in the rest of each of the above matrices (N_y_a, Seine Catch age-comp and Total Seine catch) in 1 loop... Thank you, Hulson!
      //Below fills from 1985 to 1992
      N_y_a(i,1)=init_age_0(i);
      for(int j=2;j<=nage-1;j++){
        if(turn_on_effort){
          gc(i-1,j-1)=value(gc_V(i-1,j-1)*N_y_a(i-1,j-1));
          pc(i-1,j-1)=value(pc_V(i-1,j-1)*N_y_a(i-1,j-1));
          fbc(i-1,j-1)=value(fbc_V(i-1,j-1)*(N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc(i-1,j-1)+pk*pc(i-1,j-1)))*Sur_summer(i-1,j-1));
        }

        // Still no disease data, so no dissease survival component (11/15/2020)
        N_y_a(i,j) = ((N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc(i-1,j-1)+pk*pc(i-1,j-1)))*Sur_summer(i-1,j-1)-fbc(i-1,j-1))*Sur_winter(i,j-1);

      }
      // Plus group
      int j=nage;
      if(turn_on_effort){
          gc(i-1,j)=value(gc_V(i-1,j)*N_y_a(i-1,j));
          pc(i-1,j)=value(pc_V(i-1,j)*N_y_a(i-1,j));
          fbc(i-1,j)=value(fbc_V(i-1,j)*(N_y_a(i-1,j)-(SeAC(i-1,j)*N_se(i-1)+gc(i-1,j)+pk*pc(i-1,j)))*Sur_summer(i-1,j));
      }

      N_y_a(i,j)=((N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc(i-1,j-1)+pk*pc(i-1,j-1)))*Sur_summer(i-1,j-1)-fbc(i-1,j-1))*Sur_winter(i,j-1) + ((N_y_a(i-1,j)-(SeAC(i-1,j)*N_se(i-1)+gc(i-1,j)+pk*pc(i-1,j)))*Sur_summer(i-1,j)-fbc(i-1,j))*Sur_winter(i,j);
      
      for(int j=1;j<=nage;j++){
        dvariable pen2=0.0;
        N_y_a(i,j)=posfun(N_y_a(i,j), 0.01, pen2);
        if(N_y_a(i,j)<=not_below_this){
          penCount+=1;
        }
        f_llk+=1000*pen2;
        CC_top(i,j)=N_y_a(i,j)*Vul(j);
      } 
      CC_bot=rowsum(CC_top); 
      for(int j=1;j<=nage;j++){
        SeAC(i,j)=CC_top(i,j)/CC_bot(i);
        bio_a_a(i,j)=SeAC(i,j)*w_a_a(i,j);
      }
      SeBiomass=rowsum(bio_a_a);

      if(turn_on_effort){
        // For simulation, calculate new seine yield based on specified exploitation/harvest rate sc_F from simulation section
        sc(i)=value(sc_F(i)*sum(elem_prod(Vul(1,nage),elem_prod(N_y_a(i)(1,nage),w_a_a(i)(1,nage)))));
      }
      N_se(i)=sc(i)/SeBiomass(i);

    } else if(i>13){
        //Below fills from 1993 to nyr_tobefit - matches ADF&G Excel matrix exactly
        N_y_a(i,1)=init_age_0(i);
        for(int j=2;j<=nage-1;j++){
          if(turn_on_effort){
            gc(i-1,j-1)=value(gc_V(i-1,j-1)*N_y_a(i-1,j-1));
            pc(i-1,j-1)=value(pc_V(i-1,j-1)*N_y_a(i-1,j-1));
            fbc(i-1,j-1)=value(fbc_V(i-1,j-1)*(N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc(i-1,j-1)+pk*pc(i-1,j-1)))*Sur_summer(i-1,j-1));
          }

          if((i-1)>=vhsv_start_est){
            Sur_vhsv(i-1,j-1) = (1 - Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1)) + Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1)*pseudo_rec(i-1);

            immune(i,j) = ( immune(i-1,j-1)+Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1)*pseudo_rec(i-1) ) / ( (1 - Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1))+Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1)*pseudo_rec(i-1) );
            suscep(i,j) = 1 - immune(i,j);
          }

          N_y_a(i,j)=((N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc(i-1,j-1)+pk*pc(i-1,j-1)))*Sur_summer(i-1,j-1)-fbc(i-1,j-1))*Sur_winter(i,j-1) * Sur_vhsv(i-1,j-1);

        }
        
        // Plus group
        int j=nage;
        if(turn_on_effort){
          gc(i-1,j)=value(gc_V(i-1,j)*N_y_a(i-1,j));
          pc(i-1,j)=value(pc_V(i-1,j)*N_y_a(i-1,j));
          fbc(i-1,j)=value(fbc_V(i-1,j)*(N_y_a(i-1,j)-(SeAC(i-1,j)*N_se(i-1)+gc(i-1,j)+pk*pc(i-1,j)))*Sur_summer(i-1,j));
        }
        
        dvariable Nya_prenage = ((N_y_a(i-1,j-1)-(SeAC(i-1,j-1)*N_se(i-1)+gc(i-1,j-1)+pk*pc(i-1,j-1)))*Sur_summer(i-1,j-1)-fbc(i-1,j-1))*Sur_winter(i,j-1);
        dvariable Nya_curnage = ((N_y_a(i-1,j)-(SeAC(i-1,j)*N_se(i-1)+gc(i-1,j)+pk*pc(i-1,j)))*Sur_summer(i-1,j)-fbc(i-1,j))*Sur_winter(i,j);

        if((i-1)>=vhsv_start_est){
          Sur_vhsv(i-1,j-1) = (1 - Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1)) + Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1)*pseudo_rec(i-1);
          Sur_vhsv(i-1,j) = (1 - Vul_sympat(j)*suscep(i-1,j)*annual_inf(i-1)) + Vul_sympat(j)*suscep(i-1,j)*annual_inf(i-1)*pseudo_rec(i-1);

          dvariable immune_prenage = ( immune(i-1,j-1)+Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1)*pseudo_rec(i-1) ) / ( (1 - Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1))+Vul_sympat(j-1)*suscep(i-1,j-1)*annual_inf(i-1)*pseudo_rec(i-1) );
          dvariable immune_curnage = ( immune(i-1,j)+Vul_sympat(j)*suscep(i-1,j)*annual_inf(i-1)*pseudo_rec(i-1) ) / ( (1 - Vul_sympat(j)*suscep(i-1,j)*annual_inf(i-1))+Vul_sympat(j)*suscep(i-1,j)*annual_inf(i-1)*pseudo_rec(i-1) );
          immune(i,j) = (immune_prenage * Nya_prenage + immune_curnage * Nya_curnage)/(Nya_prenage + Nya_curnage);
          suscep(i,j) = 1 - immune(i,j);
        }

        N_y_a(i,j) = Nya_prenage * Sur_vhsv(i-1,j-1) +  Nya_curnage * Sur_vhsv(i-1,j);

        // Check for impossible values (negative numbers) & blow up likelihood if present
        for(int j=1;j<=nage;j++){
          dvariable pen3=0.0;
          N_y_a(i,j)=posfun(N_y_a(i,j), 0.01, pen3);
          if(N_y_a(i,j)<=not_below_this){
            penCount+=1;
          }
          f_llk+=1000*pen3;
          CC_top(i,j)=N_y_a(i,j)*Vul(j);
        }
        CC_bot=rowsum(CC_top); 
        for(int j=1;j<=nage;j++){
          SeAC(i,j)=CC_top(i,j)/CC_bot(i);
          bio_a_a(i,j)=SeAC(i,j)*w_a_a(i,j);
        }
        SeBiomass=rowsum(bio_a_a);

        if(turn_on_effort){
          // For simulation, calculate new seine yield based on specified exploitation/harvest rate sc_F from simulation section
          sc(i)=value(sc_F(i)*sum(elem_prod(Vul(1,nage),elem_prod(N_y_a(i)(1,nage),w_a_a(i)(1,nage)))));
        }
        N_se(i)=sc(i)/SeBiomass(i);
    }

    //Pre-Fishery Spawning Biomass or Pre-Fishery Run Biomass, mt
    if(DD_Mat==1){
    	Mat(i)(1,3)=0;
	  	Mat(i,4)=1/(1+exp(matur_age3_per1+exp(matur_age4_per1)*sum(N_y_a(i)(1,nage))));
	  	Mat(i,5)=Mat(i,4)+(1-Mat(i,4))/(1+exp(matur_age3_per2));
	  	Mat(i)(6,nage)=1;
  	}
    SSB(i) = 0;

    if(mat_mod_type==1){
      for(int j=1;j<=nage;j++){
        Early_Sp_biomass(i,j)=Mat(i,j)*N_y_a(i,j)*w_a_a(i,j);
        SSB(i)+= Early_Sp_biomass(i,j); //Spawning Stock Biomass; sum over ages the pre-fishery spawning biomass by year
        //rowsum(Early_Sp_biomass);
      }
    }else{
      for(int j=1;j<=nage;j++){
        Early_Sp_biomass(i,j)=(Vul_survey(j)*Mat(i,j)+(1-Vul_survey(j))*Mat_unobs(i,j))*N_y_a(i,j)*w_a_a(i,j);
        SSB(i)+= Early_Sp_biomass(i,j); //Spawning Stock Biomass; sum over ages the pre-fishery spawning biomass by year
      }
    }

    SB_star[i]=SSB[i]*2204.62/2000; // pre-fishery run biomass in TONS

    //Pre-Fishery Spawning Age-Composition
    SpAC_bot(i) = 0;
    if(mat_mod_type==1){
      for(int j=1;j<=nage;j++){
        SpAC_top(i,j)=Mat(i,j)*N_y_a(i,j); //numerator of Spawning Age-Comp(SpAC)
        SpAC_bot(i)+=SpAC_top(i,j);
        //SpAC_bot=rowsum(SpAC_top); //denom of SpC
      }
    }else{
      for(int j=1;j<=nage;j++){
        SpAC_top(i,j)=Vul_survey(j)*N_y_a(i,j); //numerator of Spawning Age-Comp(SpAC)
        SpAC_bot(i)+=SpAC_top(i,j);
        //SpAC_bot=rowsum(SpAC_top); //denom of SpC
      }
    }

    for(int j=1;j<=nage;j++){
      SpAC(i,j)=SpAC_top(i,j)/SpAC_bot(i);  //fills in Spawning Age-Comp
    }

    // Post-First half year Fisheries Spawning Population Estimates, called Naturally Spawning Pop in Excel model
    // N_sp Spawning Population (in millions) - mature abundance with spring catch and impound removed
    // Set up the remaining age-classes each year (incrementally by year beginning with age 5 in 1980) are zero
    // First set first 4 rows all to zero
    if(i<=4){
      for(int j=1;j<=nage;j++){
        N_sp(i,j)=0;
        SPB(i,j)=0;
      }
      // Now generate naturally spawning pop as correct structure for first four years
      SB(i)=0;
      for(int j=4;j<=m;j++){
        
        N_sp(i,j)=Mat(i,j)*(N_y_a(i,j)-(SeAC(i,j)*N_se(i)+gc(i,j)+pc(i,j)));
        
        dvariable pen4=0.0;
        N_sp(i,j)=posfun(N_sp(i,j), .01, pen4);
        if(N_sp(i,j)<=not_below_this){
          penCount+=1;
        }
        f_llk+=1000*pen4;

        SPB(i,j)=N_sp(i,j)*w_a_a(i,j);
        //SPB(i,j)=N_sp(i,j)*w_a_a(i,j)/Mat(i,j);
        
        SB(i)+=SPB(i,j); // SB=rowsum(SPB); //Total naturally spawning biomass
      } 
      m++;
    } else if(i>=5){
      // Now fill in the remaining rows regularly
      SB(i)=0;
      for(int j=4;j<=nage;j++){

        N_sp(i,j)=Mat(i,j)*(N_y_a(i,j)-(SeAC(i,j)*N_se(i)+gc(i,j)+pc(i,j)));

        dvariable pen5=0.0;
        N_sp(i,j)=posfun(N_sp(i,j), .01, pen5);
        if(N_sp(i,j)<=not_below_this){
          penCount+=1;
        }
        f_llk+=1000*pen5;
        
        SPB(i,j)=N_sp(i,j)*w_a_a(i,j);
        //SPB(i,j)=N_sp(i,j)*w_a_a(i,j)/Mat(i,j);
        
        SB(i)+=SPB(i,j);
      }
    }
    
  }

  // Infection incidence, fatality, and immunity of spawning population
  inf_inc_sp.initialize();
  fatal_sp.initialize();
  seroprev_sp.initialize();

  inf_inc_age3.initialize();
  fatal_age3.initialize();
  seroprev_age3.initialize();

  dvar_vector temp_inc(1,nage);
  dvar_vector temp_fat(1,nage);
  dvar_vector temp_imm(1,nage);

  for(int i=vhsv_start_est;i<=nyr_tobefit;i++){
    for(int j=1;j<=nage;j++){
      temp_inc(j) = Mat(i,j)*N_y_a(i,j)*Vul_sympat(j)*suscep(i,j)*annual_inf(i);
      temp_fat(j) = Mat(i,j)*N_y_a(i,j)*Vul_sympat(j)*suscep(i,j)*annual_inf(i)*(1-pseudo_rec(i));
      temp_imm(j) = Mat(i,j)*N_y_a(i,j)*immune(i,j);
    }
    inf_inc_sp(i) = sum(temp_inc)/sum(N_sp(i)(1,nage));
    fatal_sp(i) = sum(temp_fat)/sum(N_sp(i)(1,nage));
    seroprev_sp(i) = sum(temp_imm)/sum(N_sp(i)(1,nage));

    inf_inc_age3(i) = (Mat(i,4)*N_y_a(i,4)*Vul_sympat(4)*suscep(i,4)*annual_inf(i))/N_sp(i,4);
    fatal_age3(i) = (Mat(i,4)*N_y_a(i,4)*Vul_sympat(4)*suscep(i,4)*annual_inf(i)*(1-pseudo_rec(i)))/N_sp(i,4);
    seroprev_age3(i) = (Mat(i,4)*N_y_a(i,4)*immune(i,4))/N_sp(i,4);
  }


FUNCTION void calc_surveyvalues()
  MDM.initialize(); 
  EggAC.initialize();
  EggAC_sum.initialize();
  EGG.initialize();
  Early_bio_a_a.initialize();
  Early_biomass.initialize();  
  HYD_ADFG.initialize();
  HYD_PWSSC.initialize();
  juv_pred.initialize();
  sero_pred.initialize();

  //Egg deposition - this data set patchily exists for 10 out of the nyr_tobefit years
  for(int i=1;i<=nyr_tobefit;i++){
    if(egg(i)==-9){
      for(int j=1;j<=nage;j++){
        EggAC(i,j)=0;
      } 
    } else{
      for(int j=1;j<=nage;j++){
        EggAC(i,j)=N_sp(i,j)*fecun(i,j);
      }
    }
  }
    
  EggAC_sum=rowsum(EggAC);
  for(int i=1;i<=nyr_tobefit;i++){
    EGG(i)=0.000001*f_sp(i)*EggAC_sum(i); //f_sp is female spawners
  }

  //Mile-days of milt
  for(int i=1;i<=nyr_tobefit;i++){
    MDM(i)=(1-f_sp(i))*SB(i)/mdm_c;
  }

  // ADFG & PWSSC Hydroacoustic Survey Biomass 
  for(int i=1;i<=nyr_tobefit;i++){
  	for(int j=1;j<=nage;j++){
      if(mat_mod_type==1){
        //Early_bio_a_a(i,j)=N_y_a(i,j)*w_a_a(i,j);
        Early_bio_a_a(i,j)=Mat(i,j)*N_y_a(i,j)*w_a_a(i,j);
      }else{
        Early_bio_a_a(i,j)=(Vul_survey(j)*Mat(i,j)+(1-Vul_survey(j))*Mat_unobs(i,j))*N_y_a(i,j)*w_a_a(i,j);
      }
    }
 	  Early_biomass=rowsum(Early_bio_a_a); // includes mature and non-mature fish
    HYD_ADFG(i)=Early_biomass(i)*exp(hydADFG_q);
    HYD_PWSSC(i)=Early_biomass(i)*exp(hydPWSSC_q);
  }

  // Aerial juvenile survey (incorporated 12/2019)
  for(int i=1;i<=nyr_tobefit;i++){
    juv_pred(i)=N_y_a(i,2)*exp(log_q_juv); 
  }

  // Seroprevalence survey values (11/2020)
  Vul_sero_survey.initialize();
  for (int a=1;a<=nage;a++){
    // Vul_sero_survey(a) = 1/(1+exp(-log(19)*((a-1)-sero_samp_vul_a50)/(sero_samp_vul_a95-sero_samp_vul_a50)));
    Vul_sero_survey(a) = 1/(1+exp(-log(19)*(3-matur_age3_per1)/(matur_age4_per1-matur_age3_per1)));
  }

  int k = 1;

  for(int i=vhsv_start_est;i<=nyr_tobefit;i++){
    for(int j=1;j<=nage;j++){
      sero_pred(i,k) = Vul_sero_survey(j)*N_y_a(i,j)*immune(i,j); // This assumes seroprevalence samples target all mature individuals (i.e. maturity = survey vulnerability)
      k+=1;
      sero_pred(i,k) = Vul_sero_survey(j)*N_y_a(i,j)*(1-immune(i,j));
      k+=1;
    } 
    k=1;
    sero_pred(i)(1,n_sero) = sero_pred(i)(1,n_sero)/sum(sero_pred(i)(1,n_sero));
  }


FUNCTION void calc_nll_components()
  age0_devs_penllk.initialize();
  mort_devs_penllk.initialize();
  age0_covar_prior.initialize();
  mort_covar_prior.initialize();

  // Remove penalized lik for unconstrained rec devs - 12/22/2019
  // if(sigma_age0devs==0){
  //   age0_devs_penllk = 0;
  // }else{
  //  for(int i=1; i<=nyr_tobefit; i++){
  //    if(R_cov_model==1){
  //      age0_devs_penllk += log(sigma_age0devs)+0.5*square(annual_age0devs(1,i))/square(sigma_age0devs);
  //    }else if(R_cov_model==2){
  //      for (int k=1;k<=rec_cov_counter;k++){
  //        if(beta_recruit_ind(k)==0){
  //        }else if(age0_covariates(i,beta_recruit_ind(k))!=-9){
  //            // age0_devs_penllk += log(sigma_age0devs)+0.5*square(annual_age0devs(i))/square(sigma_age0devs);
  //          age0_covar_prior += log(sigma_age0covar(k))+0.5*square(annual_age0devs(k,i)-age0_covariates(i,beta_recruit_ind(k)))/square(sigma_age0covar(k));
  //        }
  //      }
  //      age0_devs_penllk += log(sigma_age0devs)+0.5*square(colsum(annual_age0devs)(i))/square(sigma_age0devs);
  //    }
  //  }
  // }

  if(M_cov_model==2){
    for(int i=1; i<=nyr_tobefit; i++){
      for (int k=1;k<=mor_cov_counter;k++){
        if(beta_mortality_ind(k)==0){
        }else if(mor_covariates(i,beta_mortality_ind(k))!=-9){
          mort_covar_prior += log(sigma_morcovar(k))+0.5*square(annual_mortdevs(k,i)-mor_covariates(i,beta_mortality_ind(k)))/square(sigma_morcovar(k));
        }
      }
      mort_devs_penllk += log(sigma_mortdevs)+0.5*square(colsum(annual_mortdevs)(i))/square(sigma_mortdevs);
    }
    //cout << annual_mortdevs << endl << endl;
    //cout << colsum(annual_mortdevs) << endl; 
  }


  Setemp_1.initialize();
  Sptemp_1.initialize();
  MDMllk_ind.initialize();
  EGGllk_ind.initialize();
  H_ADFGllk_ind.initialize();
  H_PWSSCllk_ind.initialize();
  juv_llk_ind.initialize();
  sero_llk.initialize();


  //Seine Age Composition - this data set is very patchy
    for(int i=1;i<=nyr_tobefit;i++){
      for(int j=1;j<=nage;j++){
        if(seine(i,j)<=0){
          Setemp_1(i,j)=0;
        }else if(SeAC(i,j)==0){
          Setemp_1(i,j)=0;
        }else{
          Setemp_1(i,j)=seine(i,j)*(log(SeAC(i,j))-log(seine(i,j)));
        }
      }
    }
    
    Setemp_2=rowsum(Setemp_1);
    for(int i=1;i<=nyr_tobefit;i++){
      Setemp_3(i)=ESS_Se(i)*Setemp_2(i);
    }
    Se_llk=-sum(Setemp_3);
  
  //Spawning Age Composition
    for(int i=1;i<=nyr_tobefit;i++){
      for(int j=1;j<=nage;j++){
        if(spac(i,j)<=0){
          Sptemp_1(i,j)=0;
        }else if(SpAC(i,j)==0){
          Sptemp_1(i,j)=0;
        }else{
          Sptemp_1(i,j)=spac(i,j)*(log(SpAC(i,j))-log(spac(i,j)));
        }
      }
    }

    Sptemp_2=rowsum(Sptemp_1);
    for(int i=1;i<=nyr_tobefit;i++){
      Sptemp_3(i)=ESS_Sp(i)*Sptemp_2(i);
    }
    Sp_llk=-sum(Sptemp_3);

  //Mile-days of milt Likelihood component
  for(int i=1;i<=nyr_tobefit;i++) {
    MDMtemp_1(i)=log(MDM(i))-log(mdm(i));
    MDMllk_ind(i)=(MDMtemp_1(i)*MDMtemp_1(i))/(2*m_add*m_add)+log(m_add);
  }
  MDMtemp_2=norm2(MDMtemp_1);
  //M_VAR=MDMtemp_2/nyr_tobefit+(m_add*m_add);
  MDMllk=nyr_tobefit*log(m_add)+(.5*MDMtemp_2/(m_add*m_add));

  //Egg Deposition Likelihood component
  EGGllk=0;
  for(int i=1;i<=nyr_tobefit;i++){
    Eg_SD(i)=0;
    EGGtemp(i)=0;
    EGGllk_ind(i)=0;
  }
  for(int i=5;i<=5;i++){
    Eg_SD(i)=sqrt((cv_egg(i)*cv_egg(i))+(egg_add*egg_add));
    EGGtemp(i)=(log(EGG(i))-log(egg(i)));
    EGGllk_ind(i)=log(Eg_SD(i))+(.5*EGGtemp(i)*EGGtemp(i)/(Eg_SD(i)*Eg_SD(i)));
    EGGllk+=EGGllk_ind(i);
  }
  for(int i=9;i<=13;i++){
    Eg_SD(i)=sqrt((cv_egg(i)*cv_egg(i))+(egg_add*egg_add));
    EGGtemp(i)=(log(EGG(i))-log(egg(i)));
    EGGllk_ind(i)=log(Eg_SD(i))+(.5*EGGtemp(i)*EGGtemp(i)/(Eg_SD(i)*Eg_SD(i)));
    EGGllk+=EGGllk_ind(i);
  }
  for(int i=15;i<=18;i++){
    Eg_SD(i)=sqrt((cv_egg(i)*cv_egg(i))+(egg_add*egg_add));
    EGGtemp(i)=(log(EGG(i))-log(egg(i)));
    EGGllk_ind(i)=log(Eg_SD(i))+(.5*EGGtemp(i)*EGGtemp(i)/(Eg_SD(i)*Eg_SD(i)));
    EGGllk+=EGGllk_ind(i);
  }
  // EGGllk_ind = EGGllk_ind*10;
  EGGllk = EGGllk;

  //ADFG Hydroacoustic Survey Biomass Likelihood component
  for(int i=1;i<=hydADFG_start-1;i++){
    HtempADFG_vec(i)=0;
  }

  int N_hydADFG=0;
  for(int i=hydADFG_start;i<=hydADFG_start+4;i++){ 
    // hyd~_start variable holds index of first year of
    // survey depending on data source
    // UPDATED 07/23/2015 to reflect 2 additional years of missing data
    HtempADFG_vec(i)=log(hydADFG(i))-log(HYD_ADFG(i));
    H_ADFGllk_ind(i)=(HtempADFG_vec(i)*HtempADFG_vec(i))/(hydADFG_add*hydADFG_add*2)+log(hydADFG_add);
    N_hydADFG+=1;
  }
  for(int i=hydADFG_start+4+1;i<=nyr_tobefit;i++){
    HtempADFG_vec(i)=0;    //missing final 5 years of ADF&G data as of 07/2015
    H_ADFGllk_ind(i)=0;
  }
  HtempADFG_num=norm2(HtempADFG_vec);
  H_ADFGllk=(N_hydADFG)*log(hydADFG_add)+(0.5*HtempADFG_num/(hydADFG_add*hydADFG_add));
  // minus 5 to account for the missing final 5 years of ADF&G data

  //PWSSC Hydroacoustic Survey Biomass Likelihood component
  H_PWSSCllk=0;

  int N_hydPWWSC=0;
  for(int i=1;i<=nyr_tobefit;i++){ // hyd_start variable holds index of first year of survey depending on data source
    if(cv_hydPWSSC(i)==-9) {
      PWSSC_SD(i)=0;
      HtempPWSSC_vec(i)=0;
      H_PWSSCllk_ind(i)=0;
    }
    else{
      PWSSC_SD(i)=sqrt((cv_hydPWSSC(i)*cv_hydPWSSC(i))+(hydPWSSC_add*hydPWSSC_add));
      HtempPWSSC_vec(i)=(log(hydPWSSC(i))-log(HYD_PWSSC(i)));
      H_PWSSCllk_ind(i)=log(PWSSC_SD(i))+(.5*HtempPWSSC_vec(i)*HtempPWSSC_vec(i)/(PWSSC_SD(i)*PWSSC_SD(i)));
      H_PWSSCllk+=H_PWSSCllk_ind(i);
      N_hydPWWSC+=1;
    }
  }

  // Aerial juvenile survey (incorporated 12/2019)
  juv_llk=0;
  for(int i=1;i<=nyr_tobefit;i++) {
    if(juv_ind(i)!=-9) {
      juv_llk_ind(i)=dnbinom(juv_ind(i),juv_pred(i),exp(log_overdisp_juv));
      juv_llk+=juv_llk_ind(i);
    }
  }
  
  // Seroprevalence obs likelihood - multinomial (may need to come back to) (11/2020)
  
  // dvector sero_ss(1,nyr_tobefit); // This is when I had that actual sample #'s in the PWS_ASA.dat
  // sero_ss = rowsum(sero_obs);

  dvar_matrix sero_llk_ind(1,nyr_tobefit,1,n_sero);

  for(int i=sero_start;i<=nyr_tobefit;i++){
    // sero_obs_comp(i)(1,n_sero) = sero_obs(i)(1,n_sero)/sero_ss(i);

    dvar_vector sero_obs_comp(1,n_sero);
    sero_obs_comp = sero_obs(i)(1,n_sero);

    dvar_vector sero_pred_comp(1,n_sero);
    sero_pred_comp = sero_pred(i)(1,n_sero);

    dvariable x_bin;
    dvariable n_bin;
    dvariable p_bin;

    for(int j=1;j<=nage;j++){
      if(sero_obs_comp(j*2-1)<=0 | sero_pred_comp(j*2-1)==0){
        sero_llk_ind(i,j)=0;
      }else{

        // Binomial
        x_bin = sero_obs_comp(j*2-1)*ESS_Antibody(i);
        n_bin = (sero_obs_comp(j*2-1)+sero_obs_comp(j*2))*ESS_Antibody(i);
        p_bin = sero_pred_comp(j*2-1)/(sero_pred_comp(j*2-1) + sero_pred_comp(j*2));
        sero_llk_ind(i,j) = dbinom(x_bin,n_bin,p_bin);

        // Multinomial
        // sero_llk_ind(i,j) = sero_obs_comp(j)*(log(sero_pred_comp(j)) - log(sero_obs_comp(j))); // elem_prod(sero_obs_comp(i)(1,n_sero),log(sero_pred(i)(1,n_sero)) - log(sero_obs_comp(i)(1,n_sero)));
      }
    }
    // cout << sero_llk_ind << endl;
  }

  // Binomial
  sero_llk = sum(rowsum(sero_llk_ind));

  // Multinomial
  // sero_llk = -sum(elem_prod(ESS_Antibody,rowsum(sero_llk_ind)));
  // sero_llk = 0;

  Z_prior = 0;

  hydADFG_add_prior = log(0.03)+0.5*square((hydADFG_add-0.3)/(0.03));
  hydPWSSC_add_prior = log(0.03)+0.5*square((hydPWSSC_add-0.3)/(0.03));
  m_add_prior = log(0.03)+0.5*square((m_add-0.3)/(0.03));


FUNCTION project_biomass    
  //Use the average differences across the last 5 years...
  for (int a=1; a<=nage; a++){
    tempWgt(a) = 0;
    for (int i=0; i<=4; i++){
      tempWgt(a) += w_a_a(nyr_tobefit-i, a);
    }
    avgWgt5Yr(a) = tempWgt(a)/5;
  }
  
  // Now using mean of log-recruits for projection; this should be near the median if the recruits are log-Normally distributed.
  //cout << init_age_0(33)<< endl;
  temp1MeanLgRec = 1;
  for(int i=nyr_tobefit-9;i<=nyr_tobefit;i++){
    temp1MeanLgRec *= N_y_a(i,4);
  }
  temp2MeanLgRec = log(temp1MeanLgRec)/10;
  meanLgRec = exp(temp2MeanLgRec);
  //cout << endl << "meanLgRec: " << meanLgRec << " "<< endl;

  // Plug age-3 biomass into first element of N_y_a vector for troubleshooting and display in report file
  projected_N_y_a(4) = meanLgRec;
  
  // Plug age-3 biomass into first element of vector for calculations
  projected_Early_Sp_biomass(4) = Mat(nyr_tobefit,4)*meanLgRec*avgWgt5Yr(4);
  
  // Fill in half-year survival rates for forecast year
  for(int j=1;j<=nage;j++){
    forecast_Sur_winter(j)=exp(-(0.5*M_change(nyr_tobefit)*Z_0_8offset+0.5*Z_0_8+forecast_winter_effect(j)));
  }

  // Call numbers this year using last year's info for ages 4-8
  for(int j=5;j<=nage-1;j++){
    projected_N_y_a(j)=((N_y_a(nyr_tobefit,j-1)-(SeAC(nyr_tobefit,j-1)*N_se(nyr_tobefit)+gc(nyr_tobefit,j-1)+pk*pc(nyr_tobefit,j-1)))*Sur_summer(nyr_tobefit,j-1)-fbc(nyr_tobefit,j-1))*forecast_Sur_winter(j-1);
  }
  
  // Calc numbers this year using last year's info for plus group
  for(int j=nage;j<=nage;j++){
    projected_N_y_a(j)=((N_y_a(nyr_tobefit,j-1)-(SeAC(nyr_tobefit,j-1)*N_se(nyr_tobefit)+gc(nyr_tobefit,j-1)+pk*pc(nyr_tobefit,j-1)))*Sur_summer(nyr_tobefit,j-1)-fbc(nyr_tobefit,j-1))*forecast_Sur_winter(j-1)+((N_y_a(nyr_tobefit,j)-(SeAC(nyr_tobefit,j)*N_se(nyr_tobefit)+gc(nyr_tobefit,j)+pk*pc(nyr_tobefit,j)))*Sur_summer(nyr_tobefit,j)-fbc(nyr_tobefit,j))*forecast_Sur_winter(j);
  }
  
  // Make it biomass
  for(int j=5;j<=nage;j++){
    projected_Early_Sp_biomass(j) = Mat(nyr_tobefit,j)*projected_N_y_a(j)*avgWgt5Yr(j);
  }
      
  // Take total pre-fishery biomass for projection year
  projected_PFRB = sum(projected_Early_Sp_biomass);


FUNCTION write_chain_results

    struct stat buf;
    if(stat("mcmc_out",&buf)!=0){
      system("mkdir mcmc_out");
    }

    ofstream MCMCreport1("mcmc_out/VarsReport.csv",ios::app);
    ofstream MCMCreport2("mcmc_out/Age3.csv",ios::app);
    ofstream MCMCreport3("mcmc_out/HYD_ADFG.csv",ios::app);
    ofstream MCMCreport4("mcmc_out/HYD_PWSSC.csv",ios::app);
    ofstream MCMCreport5("mcmc_out/EGG.csv",ios::app);
    ofstream MCMCreport6("mcmc_out/MDM.csv",ios::app);
    ofstream MCMCreport7("mcmc_out/PostFRbiomass.csv",ios::app);
    ofstream MCMCreport8("mcmc_out/SeAC.csv",ios::app); // writes Seine age comps from each iteration to a file
    ofstream MCMCreport9("mcmc_out/SpAC.csv",ios::app); // writes spawner age comps from each iteration to a file
    ofstream MCMCreport10("mcmc_out/juv_schools.csv",ios::app); 
    ofstream MCMCreport11("mcmc_out/Num_at_age.csv",ios::app); 
    ofstream MCMCreport12("mcmc_out/Sero_pred.csv",ios::app); 
    ofstream parReport("mcmc_out/iterations.csv",ios::app);
    ofstream LLikReport("mcmc_out/llikcomponents.csv",ios::app);
    ofstream PFRReport("mcmc_out/PFRBiomass.csv",ios::app);
    ofstream indiv_LLikReport("mcmc_out/llik_observations.csv",ios::app);
    ofstream recruit_effect_report("mcmc_out/recruitment_effects.csv",ios::app);
    ofstream summer_survival_report("mcmc_out/adult_survival_effects_summer.csv",ios::app);
    ofstream winter_survival_report("mcmc_out/adult_survival_effects_winter.csv",ios::app);
    ofstream vhsv_agespec_sur_report("mcmc_out/vhsv_survival.csv",ios::app);
    ofstream vhsv_infec_incid_report("mcmc_out/vhsv_infection.csv",ios::app);
    ofstream vhsv_seroprev_report("mcmc_out/vhsv_seroprev.csv",ios::app);
    ofstream vhsv_fatal_rate_report("mcmc_out/vhsv_fatality.csv",ios::app);
    ofstream vhsv_infec_incid_age3_report("mcmc_out/vhsv_infection_age3.csv",ios::app);
    ofstream vhsv_seroprev_age3_report("mcmc_out/vhsv_seroprev_age3.csv",ios::app);
    ofstream vhsv_fatal_rate_age3_report("mcmc_out/vhsv_fatality_age3.csv",ios::app);

    MCMCreport1 << m_add <<  "," << egg_add << "," << hydADFG_add  << "," << hydPWSSC_add << endl;
  
    for (int i=1; i<=nyr_tobefit-1; i++){
        MCMCreport2 << N_y_a(i,4) << ","; 
        }
    MCMCreport2 << N_y_a(nyr_tobefit,4) << endl; // this is the projected recruitment for the latest year
    
    for (int i=1; i<=nyr_tobefit-1; i++){
        MCMCreport3 << HYD_ADFG(i) << ","; 
        }
    MCMCreport3 << HYD_ADFG(nyr_tobefit) << endl;
    
    for (int i=1; i<=nyr_tobefit-1; i++){
        MCMCreport4 << HYD_PWSSC(i) << ","; 
        }
    MCMCreport4 << HYD_PWSSC(nyr_tobefit) << endl;
    
    for (int i=1; i<=nyr_tobefit-1; i++){
        MCMCreport5 << EGG(i) << ","; 
    }
    MCMCreport5 << EGG(nyr_tobefit) << endl;
  
    for (int i=1; i<=nyr_tobefit-1; i++){
        MCMCreport6 << MDM(i) << ","; 
    }
    MCMCreport6 << MDM(nyr_tobefit) << endl;
    
    // SB is Naturally spawning biomass
    for (int i=1; i<=nyr_tobefit-1; i++){
        MCMCreport7 << SB(i) << ","; 
    } 
    MCMCreport7 << SB(nyr_tobefit) << endl;
  
    // write age comps (first Seine, then Spawner) to a file
    for (int i=1; i<=nyr_tobefit-1; i++){
      for (int j=1; j<=nage; j++){
        MCMCreport8 << SeAC(i,j) << ",";
      }
    }
    
    for (int j=1; j<=nage-1; j++){
        MCMCreport8 << SeAC(nyr_tobefit,j) << ","; 
    }
    MCMCreport8 << SeAC(nyr_tobefit,nage) << endl;
   
    for (int i=1; i<=nyr_tobefit-1; i++){
      for (int j=1; j<=nage; j++){
          MCMCreport9 << SpAC(i,j) << ",";
      }
    }
    
    for (int j=1; j<=nage-1; j++){
        MCMCreport9 << SpAC(nyr_tobefit,j) << ","; 
    }
    MCMCreport9 << SpAC(nyr_tobefit,nage) << endl;

    // Aerial juvenile survey (incorporated 12/2019)
    for (int i=1; i<=nyr_tobefit-1; i++){
      MCMCreport10 << juv_pred(i) << ","; 
      }
    MCMCreport10 << juv_pred(nyr_tobefit) << endl;


    // Estimated Numbers at age (incorporated 01/2020)
    for (int i=1; i<=nyr_tobefit-1; i++){
      for (int j=1; j<=nage; j++){
          MCMCreport11 << N_y_a(i,j) << ",";
      }
    }
    
    for (int j=1; j<=nage-1; j++){
        MCMCreport11 << N_y_a(nyr_tobefit,j) << ","; 
    }
    MCMCreport11 << N_y_a(nyr_tobefit,nage) << endl;

    for (int i=1; i<=nyr_tobefit-1; i++){
      for (int k=1; k<=n_sero; k++){
          MCMCreport12 << sero_pred(i,k) << ",";
      }
    }
    for (int k=1; k<=n_sero-1; k++){
        MCMCreport12 << sero_pred(nyr_tobefit,k) << ","; 
    }
    MCMCreport12 << sero_pred(nyr_tobefit,n_sero) << endl;


    // parReport records all param values for each saved iteration to a .csv
      parReport << VHSV_age3_4_mort_93 << "," << ICH_age5_8_mort_93 << "," << Z_0_8 << "," << Z_9 << ",";
      parReport << Z_0_8offset << "," << Z_9offset << ",";
      parReport << matur_age3_per1 << "," << matur_age4_per1 << "," << matur_age3_per2 << "," << matur_age4_per2 << ",";
      parReport << alpha_v << "," << beta_v << "," << survey_vul_alpha << "," << survey_vul_beta << ",";
      parReport << loginit_pop(1) << "," << loginit_pop(2) << "," << loginit_pop(3) << "," << loginit_pop(4) << "," << loginit_pop(5) << ",";
      parReport << egg_add << "," << logmdm_c << "," << m_add << ",";
      parReport << hydADFG_q << "," << hydADFG_add << ","<< hydPWSSC_q << ","<< hydPWSSC_add << ",";

      // Aerial juvenile survey (incorporated 12/2019)
      parReport << log_q_juv << "," << log_overdisp_juv << ",";
      
      for (int j=1; j<=rec_cov_counter; j++){
        for (int i=1; i<=(nyr_tobefit-3); i++){
          parReport << annual_age0devs(j,i) << ",";
        }
      }

      parReport << log_MeanAge0 << "," << Mean_Age0offset << "," << sigma_age0devs << ",";

      for (int i=1; i<=rec_cov_counter; i++){
        //switch(age0_turn_on(i)){
        //  case 1:
            parReport << beta_age0(i) << ",";
        // break;
        //} 
      }

      for (int i=1; i<=mor_cov_counter; i++){
            parReport << beta_mortality(i) << ",";
      }

      for (int j=1; j<=mor_cov_counter; j++){
        for (int i=1; i<=nyr_tobefit; i++){
          parReport << annual_mortdevs(j,i) << ",";
        }
      }

      parReport << sigma_mortdevs << ",";

      // Any offset parameters to the dev parameters for time blocks on mortality or age 9's
        for (int i=1; i<=rec_cov_counter; i++){
              parReport << beta_age0_offset(i) << ",";
        }

        for (int i=1; i<=mor_cov_counter; i++){
              parReport << beta_mortality_offset(i) << ",";
        }

      for (int i=1; i<=rec_cov_counter_age0devs; i++){
        parReport << sigma_age0covar(i) << ",";
      }

      for (int i=1; i<=mor_cov_counter; i++){
        parReport << sigma_morcovar(i) << ",";
      }

      // VHSV seroprevalence related parameters (11/2020)
      for (int i=vhsv_start_est; i<=nyr_tobefit; i++){
        parReport << annual_inf(i) << ",";
      }

      for (int i=vhsv_start_est; i<=nyr_tobefit; i++){
        parReport << recov_prob(i) << ",";
      }

      parReport << infec_vul_a50 << "," << infec_vul_a95 << "," << sero_samp_vul_a50 << "," << sero_samp_vul_a95 << ",";

      parReport << f_llk << endl;


    // Now output the loglikelihood components
    LLikReport << -Se_llk << "," << -Sp_llk << "," << EGGllk << "," << H_ADFGllk << "," << H_PWSSCllk << "," << MDMllk << "," << age0_devs_penllk << "," << mort_devs_penllk << ",";
    LLikReport << age0_covar_prior << "," << mort_covar_prior << "," << Z_prior << "," << hydADFG_add_prior << "," << hydPWSSC_add_prior << "," << m_add_prior << ",";

    // Aerial juvenile survey (12/2019) and seroprevalence samples (11/2020)
    LLikReport << juv_llk << "," << sero_llk << "," << f_llk << endl;

    indiv_LLikReport << age0_devs_penllk << "," << mort_devs_penllk << "," << Z_prior << "," << hydADFG_add_prior << "," << hydPWSSC_add_prior << "," << m_add_prior << ",";
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << -Setemp_3(i) << ",";
    }
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << -Sptemp_3(i) << ",";
    }
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << MDMllk_ind(i) << ",";
    }
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << EGGllk_ind(i) << ",";
    }
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << H_ADFGllk_ind(i) << ",";
    }
    for (int i=1; i<=(nyr_tobefit); i++){
        indiv_LLikReport << H_PWSSCllk_ind(i) << ",";
    }
    for (int i=1; i<=(nyr_tobefit-1); i++){
        indiv_LLikReport << juv_llk_ind(i) << ",";
    }
    indiv_LLikReport << juv_llk_ind(nyr_tobefit) << endl;

    for (int i=1; i<=(nyr_tobefit-1); i++){
        recruit_effect_report << exp(age0_effect(i))*Mean_Age0(i) << ",";
    }
    recruit_effect_report << exp(age0_effect(nyr_tobefit))*Mean_Age0(nyr_tobefit) << endl;
    
    for (int i=1; i<=nyr_tobefit-1; i++){
      for (int j=1; j<=nage; j++){
          summer_survival_report << Sur_summer(i,j) << ",";
      }
    }
    for (int j=1; j<=nage-1; j++){
        summer_survival_report << Sur_summer(nyr_tobefit,j) << ","; 
    }
    summer_survival_report << Sur_summer(nyr_tobefit,nage) << endl;

    for (int i=1; i<=nyr_tobefit-1; i++){
      for (int j=1; j<=nage; j++){
          winter_survival_report << Sur_winter(i,j) << ",";
      }
    }
    for (int j=1; j<=nage-1; j++){
        winter_survival_report << Sur_winter(nyr_tobefit,j) << ","; 
    }
    winter_survival_report << Sur_winter(nyr_tobefit,nage) << endl;
    

    for (int i=1; i<=nyr_tobefit-1; i++){
      for (int j=1; j<=nage; j++){
          vhsv_agespec_sur_report << Sur_vhsv(i,j) << ",";
      }
    }
    for (int j=1; j<=nage-1; j++){
        vhsv_agespec_sur_report << Sur_vhsv(nyr_tobefit,j) << ","; 
    }
    vhsv_agespec_sur_report << Sur_vhsv(nyr_tobefit,nage) << endl;


    // Disease rates on all spawners (numbers)
    for (int i=1; i<=(nyr_tobefit-1); i++){
      vhsv_infec_incid_report << inf_inc_sp(i) << ",";
    }
    vhsv_infec_incid_report << inf_inc_sp(nyr_tobefit) << endl;

    for (int i=1; i<=(nyr_tobefit-1); i++){
      vhsv_seroprev_report << seroprev_sp(i) << ",";
    }
    vhsv_seroprev_report << seroprev_sp(nyr_tobefit) << endl;

    for (int i=1; i<=(nyr_tobefit-1); i++){
      vhsv_fatal_rate_report << fatal_sp(i) << ",";
    }
    vhsv_fatal_rate_report << fatal_sp(nyr_tobefit) << endl;


    // Now disease rates on just age 3 numbers
    for (int i=1; i<=(nyr_tobefit-1); i++){
      vhsv_infec_incid_age3_report << inf_inc_age3(i) << ",";
    }
    vhsv_infec_incid_age3_report << inf_inc_age3(nyr_tobefit) << endl;

    for (int i=1; i<=(nyr_tobefit-1); i++){
      vhsv_seroprev_age3_report << seroprev_age3(i) << ",";
    }
    vhsv_seroprev_age3_report << seroprev_age3(nyr_tobefit) << endl;

    for (int i=1; i<=(nyr_tobefit-1); i++){
      vhsv_fatal_rate_age3_report << fatal_age3(i) << ",";
    }
    vhsv_fatal_rate_age3_report << fatal_age3(nyr_tobefit) << endl;
    

    // Output PFRunBiomass for each saved iteration to .csv
    for (int i=1; i<=nyr_tobefit-1; i++){
       PFRReport << SSB(i) << ","; 
       }
       PFRReport << SSB(nyr_tobefit)<< ",";
       PFRReport << projected_PFRB << endl; // Projected pre-fishery run biomass for the upcoming year

    if(DD_Mat==1){
    	ofstream maturity_report("density_dependent_maturity.csv",ios::app);
    	for (int i=1; i<=nyr_tobefit; i++){
	      for (int j=1; j<=nage; j++){
	        maturity_report << Mat(i,j) << ",";
	      }
    	}
	}

REPORT_SECTION
  // report << "foo= " << endl << setprecision(10) << foo << endl; // to set precision just for foo OR
  // report.precision(10); // in first line of report section and to set precision of all report output
  // In order to get labels at the top of each of the following .csv I open them down here, write it once, 
  // then use ::app inside the mceval phase ios

  // User defined files
  ofstream SeACreport("SeAC_pd.rep",ios::trunc);
  SeACreport << SeAC << endl;
  SeACreport.close();
  
  ofstream SpACreport("SpAC_pd.rep",ios::trunc);
  SpACreport << SpAC << endl;
  SpACreport.close();

  ofstream SeroACreport("SeroAC_pd.rep",ios::trunc);
  SeroACreport << sero_pred << endl;
  SeroACreport.close();

  // the output below is sloppy - I need to disply output for pre-1992 mort vectors without the disease index since that part of the vectors are place holders
  report<<"LOG-LIKELIHOOD COMPONENTS" << endl;
  report<< "penalty Count" << endl << penCount<< endl;
  report<<"LL OF" << endl << f_llk << endl;
  report<<"Se_llk" << endl << Se_llk <<endl;
  report<<"Sp_llk" << endl << Sp_llk <<endl;
  report<<"EGGllk " << endl << EGGllk <<endl;
  report<<"H_ADFGllk " << endl << H_ADFGllk <<endl;
  report<<"H_PWSSCllk " << endl << H_PWSSCllk <<endl;
  report<<"MDMllk " << endl <<MDMllk<<endl;
  report<<"age0_devs_penllk " << endl << age0_devs_penllk <<endl;
  report<<"mort_devs_penllk " << endl << mort_devs_penllk <<endl;
  report<<"age0_covar_prior " << endl << age0_covar_prior <<endl;
  report<<"mort_covar_prior " << endl << mort_covar_prior <<endl;
  report<<"Z_prior " << endl << Z_prior <<endl;
  report<<"hydADFG_add_prior " << endl << hydADFG_add_prior <<endl;
  report<<"hydPWSSC_add_prior " << endl << hydPWSSC_add_prior <<endl;
  report<<"m_add_prior " << endl << m_add_prior << endl;

  // Aerial juvenile survey (incorporated 12/2019)
  report<<"juv_llk " << endl << juv_llk <<endl;

  // VHSV seroprevalence (11/2020)
  report<<"sero_llk " << endl << sero_llk <<endl << endl;


  report<<"RESIDUALS" << endl;
  report<<"Seine comps residuals" << endl << Setemp_1 << endl;
  report<<"Spawner comps residuals" << endl << Sptemp_1 << endl;
  report<<"Mile-days milt residuals" << endl << MDMtemp_1 << endl;
  report<<"Egg deposition residuals" << endl << EGGtemp << endl;
  report<<"ADFG Hydroacoustic residuals" << endl << HtempADFG_vec << endl;
  report<<"PWSSC Hydroacoustic residuals" << endl << HtempPWSSC_vec << endl << endl;

  report<<"ANALYTICAL SIGMAS" << endl;
  report<<"Combined Egg SD (Eg_SD)" << endl <<Eg_SD <<endl;
  report<<"(Annual seine residuals)X(ESS)" << endl << Setemp_3 <<endl;
  report<<"(Annual spawner residuals)X(ESS)" << endl << Sptemp_3 << endl << endl;

  
  report << "DERIVED QUANTITIES" << endl;
  report << "Pre-Fishery Run Biomass in mt" << endl << SSB << endl;
  report << "Pre-Fishery Run Biomass in TONS" << endl << SB_star << endl;
  report << "Post-Fishery Spawning Biomass" << endl << SB << endl;
  report << "Estimated ADFG Hydro-acoustic Biomass" << endl << HYD_ADFG << endl;
  report << "Estimated PWSSC Hydro-acoustic Biomass" << endl << HYD_PWSSC << endl;
  report << "Pre-fishery total abundance (N_y_a)" << endl << N_y_a << endl;
  report << "Number of spawners (N_sp)" << endl << N_sp << endl << endl;

  report << "RECRUITMENT" << endl;
  report << "Recruits age-3" << endl;
  for (int i=1; i<=nyr_tobefit; i++){
    report << N_y_a(i,4) << endl; 
  }
  report << endl;

  report << "MATURITY" << endl;
  report << "Maturity-at-age of observed schools" << endl << Mat << endl << endl;
  report << "Maturity-at-age of unobserved schools" << endl << Mat_unobs << endl << endl;

  report << "ADULT SURVIVAL OUTPUTS" << endl;
  report << "Adult summer survival (Sur_summer)" << endl << Sur_summer << endl;
  report << "Adult winter survival (Sur_winter)" << endl << Sur_winter << endl << endl;
 
  report << "COVARIATE EFFECTS ON SURVIVAL" << endl;
  report << "Summer mortatlity" << endl << summer_effect << endl;
  report << "Winter mortality" << endl << winter_effect << endl << endl;

  report << "SUMMED ANNUAL MORTALITY DEVIATES (NON-ZERO IF ESTIMATED)" << endl;
  report << colsum(annual_mortdevs) << endl << endl;

  report << "ANNUAL MORTALITY DEVIATES (NON-ZERO IF ESTIMATED) - A MATRIX" << endl;
  report << annual_mortdevs << endl << endl;

  report << "VHSV ASSOCIATED SURVIVAL" << endl;
  report << "Proportion that avoids or survives infection:" << endl << Sur_vhsv << endl << endl;

  report << "VHSV AGE-SPECIFIC IMMUNITY" << endl;
  report << "Proportions that are immune to reinfection:" << endl << immune << endl << endl;

  report << "VHSV OUTBREAK CHARACTERISTICS" << endl;
  report << "Incidence rate of infection in spawning population (numbers):" << endl << inf_inc_sp << endl << endl;
  report << "Fatality rate of infection in spawning population (numbers):" << endl << fatal_sp << endl << endl;
  report << "VHSV seroprevalence in spawning population (numbers):" << endl << seroprev_sp << endl << endl;

  //report << "PROJECTED MANAGEMENT QUANTITIES" << endl;
  //report << "Mean recruits from past 10 years" << endl << meanLgRec << endl;
  //report << "Projected total pre-fishery biomass" << endl << projected_PFRB << endl;
  //report << "Projected pre-fishery biomass by age" << endl << projected_Early_Sp_biomass << endl;
  //report << "Projected numbers at age" << endl << projected_N_y_a << endl << endl;

RUNTIME_SECTION
maximum_function_evaluations 1000,10000
