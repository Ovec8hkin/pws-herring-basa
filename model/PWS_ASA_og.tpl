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
//                                                                            //
//                             Trevor A. Branch                               //
//                              tbranch@uw.edu                                //
//                                                                            //
//                               Joshua Zahner                                //
//                              jzahner @uw.edu                               //
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
    #include <statsLib.h>

    // Following adapted from thread on ADMB users Google Group
    // https://groups.google.com/forum/#!topic/admb-users/WcrSmZc_igw
    dvector rmultinom(const int& seed, const int& size, const dvector& prob)
    {  
        //Returns a multinomial sample, of size n, based on sampling probabilities p.
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
        for(i=lb+1;i<=ub;i++){
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
    // |---------------------------------------------------------------------------|
    int rseed;
    int no_estimation;
    int pin_write;
    int clean;

 LOCAL_CALCS
    adstring input_dir=".";     // directory containing input files. Defaults to the same directory as the .TPL file.
    int on = 0;
    rseed  = 0;
    no_estimation = 0;
    clean = 0;                  // Flag to delete intermediate files from output directory after run is complete.

    if (ad_comm::argc > 1){
        int on = 0;

        if((on=option_match(ad_comm::argc,ad_comm::argv,"-noest")) > -1){
            no_estimation = 1;
        }

        if((on=option_match(ad_comm::argc,ad_comm::argv,"-pinwrite")) > -1){
            pin_write = 1;
        }

        if ( (on=option_match(ad_comm::argc,ad_comm::argv,"--input_dir")) > -1 ){
            input_dir = argv[on+1];
        }

        if ( (on=option_match(ad_comm::argc,ad_comm::argv,"--clean")) > -1 ){
            clean = 1;
        }
    }

 END_CALCS

    // |---------------------------------------------------------------------------|
    // | Read in data inputs
    // |---------------------------------------------------------------------------|
    !! adstring dat_file = input_dir + "/PWS_ASA.dat";
    !! ad_comm::change_datafile_name(dat_file); 

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
    init_matrix weight_at_age(1,nyr,1,nage)
    
    //Fecundity-at-age
    init_matrix fecundity(1,nyr,1,nage)
    
    //Pound Catch
    init_matrix pound_catch(1,nyr,1,nage)
    
    //Proportion of pound fish killed
    init_number pk
    
    //Food/Bait Catch
    init_matrix food_bait_catch(1,nyr,1,nage)
    
    //Gillnet catch
    init_matrix gillnet_catch(1,nyr,1,nage)
    
    //Total seine yield
    init_vector seine_yield(1,nyr)
    
    //% Female Spawners
    init_vector female_spawners(1,nyr)
    
    //MDM
    init_vector mile_days_milt(1,nyr)
        
    //Egg Deposition
    init_vector egg_dep(1,nyr)
        
    //Standard errors derived from confidence intervals given for Egg Deposition ln(s.e.)
    init_vector cv_egg(1,nyr)       
        
    //ADFG Hydroacoustic Survey data - is a combination of ADFG and PWSSC estimates until 1994, from '95 onward are only ADFG survey data
    init_number hydADFG_start //read in year from the data file as first year of the survey: 1995
    init_vector ADFG_hydro(1,nyr)
        
    //PWSSC Hydroacoustic Survey data 
    init_number hydPWSSC_start //read in year  from the data file as first year of the survey: 1993
    init_vector PWSSC_hydro(1,nyr)
        
    //Standard errors derived from confidence intervals given for PWSSC hydroacoustic biomass ln(s.e.)
    init_vector PWSSC_hydro_cv(1,nyr)

    //Seine age distribution
    init_matrix purse_seine_age_comp(1,nyr,1,nage)

    //Spawning age composition
    init_matrix spawner_age_comp(1,nyr,1,nage)

    // Aerial juvenile survey (incorporated 12/2019)
    init_vector juvenile_survey_index(1,nyr)

    // |---------------------------------------------------------------------------|
    // | Read in VHSV seroprevalence & I. hoferi infection prevalence
    // |---------------------------------------------------------------------------|
    !! adstring disease_dat_file = input_dir + "/PWS_ASA_disease.dat";
    !! ad_comm::change_datafile_name(disease_dat_file);  

    // Seroprevalence: nyr x 2*nage matrix with each column representing # of positive then negative samples in each age (e.g. # positive age 0, # negative age 0, etc.)
    int n_age2
    !! n_age2=2*nage;
    init_matrix vhs_obs(1,nyr,1,n_age2)  
    init_number vhs_start           //start year of observations for fitting
    init_number vhs_start_est       //start year for estimating infection (can be earlier than vhs_start)
    init_int vhs_rec_cotv	        // Recovery Constant Or Time-Varing (cotv)

    // Ichthy. infection: nyr x 2*nage matrix with each column representing # of positive then negative samples in each age (e.g. # positive age 0, # negative age 0, etc.)
    init_matrix ich_obs(1,nyr,1,n_age2)  
    init_number ich_start           //start year of observations for fitting
    init_number ich_start_est       //start year for estimating infection (can be earlier than vhs_start)
    init_int ich_rec_cotv

    // |---------------------------------------------------------------------------|
    // | Read in effective sample sizes for age-composition (calculated reiteratively before running MCMC)
    // |---------------------------------------------------------------------------|
    !! adstring ess_file = input_dir + "/PWS_ASA(ESS).ctl";
    !! ad_comm::change_datafile_name(ess_file);
    init_int ESS_est
    
    !! if (ESS_est == 1) {
    !!    adstring ess_file = input_dir + "/PWS_ASA(ESS_estimate).ctl";
    !!    ad_comm::change_datafile_name(ess_file);
    !! }

    init_vector seine_ess(1,nyr_tobefit)            // Seine Effective Sample Size
    init_vector spawner_ess(1,nyr_tobefit)          // Spawing Effective Sample Size
    init_vector vhsv_antibody_ess(1,nyr_tobefit)    // Antibody sample size
    init_vector ich_ess(1,nyr_tobefit)              // Ichthyophonus sample size

    // |---------------------------------------------------------------------------|
    // | Read in the recruitment and natural mortality deviate information
    // |---------------------------------------------------------------------------|
    !! adstring covariate_file = input_dir + "/PWS_ASA(covariate).ctl";
    !! ad_comm::change_datafile_name(covariate_file);   

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
    // | Read in parameter ctl (inits,bounds,phases,priors)
    // |---------------------------------------------------------------------------|
    !! adstring par_file = input_dir + "/PWS_ASA(par).ctl";
    !! ad_comm::change_datafile_name(par_file);

    // Single value parameters
    init_int npar
    init_matrix pars_1(1,npar,1,8)

    // Vector and matrix parameers
    init_int npar2
    int totpar
    int npar3
    !! npar3 = 4*npar2; totpar = npar + npar2;
    init_ivector npar2_rag(1,npar3)
    init_matrix pars_2(1,npar3,1,npar2_rag)  // This is a ragged array
    // !! cout << pars_2 << endl;

    matrix pars_2_ctl(1,npar2,2,8)
    vector pars_2_dim(1,npar2)
    vector pars_temp(1,7)
    ivector minind2(1,npar2)
    ivector maxind2(1,npar2)
    ivector minind3(1,npar2)
    ivector maxind3(1,npar2)
    ivector j(1,4)

    int mortality_covariate_counter
    int mortality_covariate_model

    int recruitment_covariate_counter
    int recruitment_covariate_model
    int rec_cov_counter_age0devs

    int nyr_recdevs
    int nyr_vhs
    int nyr_ich

    // Maturity model - selects which one to use in the ASA (see code below)
    // 1 estimates proportions 3 & 4 mature directly from 1980-present
    // 2 estimates proportions 3 & 4 mature directly for 2 periods
    // 3 estimates maturity as a logistic function
    // 4 estimates spawner survey selectivity and fixes maturity proportions
    int maturity_model_type
    !! maturity_model_type = pars_1(5,8);


 LOCAL_CALCS
    // Design matrix for these parameters are setup similarly to pars_1 (except missing 1st column)
    // Fill in ivectors used to specify indices for ragged 3D matrix of initial values
    j.fill_seqadd(1,1);
    int i=1;
    for(int i=1; i<=npar2; i++){
        pars_temp = pars_2(j(3));
        ++pars_temp;
        pars_2_ctl(i) = pars_temp;
        pars_2_dim(i) = pars_2(j(1),1);
        if(pars_2_dim(i)==1){
            minind2(i) = pars_2(j(2),1);
            maxind2(i) = pars_2(j(2),2);
            minind3(i) = 1;
            maxind3(i) = 1;
        }else if(pars_2_dim(i)==2){
            minind2(i) = pars_2(j(2),1);
            maxind2(i) = pars_2(j(2),2);
            minind3(i) = pars_2(j(2),3);
            maxind3(i) = pars_2(j(2),4);
        }
        --pars_temp;
        j+=4;
    }
    // cout << pars_2_ctl << endl;


    // This sets up the estimable recruit and mortality covariate parameters
    // If mortality covariates are not included, counter is set to one so program does not crash (used for specifying size of estimable vector)
    // If mortality covariates are included, beta and/or deviate parameters on mortality are turned ON
    // M_cov_model default is 1, whether or not covariates are included in the model
    // M_cov_model changes to 2 if even just one of the covariates is being modeled as
    // an index that is a prior in the model
    mortality_covariate_model=1;
    for(int i=1; i<=n_mor_covs; i++){
        if(M_fixed_or_est(i)*mor_turn_on(i)==2){
            mortality_covariate_model=2;
        }
    }

    mortality_covariate_counter=sum(mor_turn_on);
    if(mortality_covariate_counter==0){
        mortality_covariate_counter=1;
    }else{
        pars_2_ctl(4,4)=2;
        //ph_betamortality=2;
        if(mortality_covariate_model==2){
            pars_2_ctl(5,4)=2;
            //ph_annualmortdevs=2;
        }
    }

    // recruitment_covariate_model default is 1, whether or not covariates are included in the model
    // recruitment_covariate_model changes to 2 if even just one of the covariates is being modeled as
    // an index that is a prior in the model
    recruitment_covariate_model=1;
    for(int i=1; i<=n_age0_covs; i++){
        if(R_fixed_or_est(i)*age0_turn_on(i)==2){
            recruitment_covariate_model=2;
        }
    }

    recruitment_covariate_counter=sum(age0_turn_on);
    if(recruitment_covariate_counter==0){
        recruitment_covariate_counter=1;
    }else{
        pars_2_ctl(3,4)=2;
        // ph_betaage0=2;
    }

    if(recruitment_covariate_model==1){
        rec_cov_counter_age0devs=1;
    }else if(recruitment_covariate_model==2){
        rec_cov_counter_age0devs=recruitment_covariate_counter;
    }


    // Adjust maxind for parameters that are determined internally here (e.g. by nyr_tobeit,recruitment_covariate_counter, etc.)
    // This probably could be improved...
    nyr_recdevs = nyr_tobefit-1;
    if(juvenile_survey_index(nyr_tobefit) == -9 || juvenile_survey_index(nyr_tobefit-1) == -9){
        nyr_recdevs = nyr_tobefit-3;
    }
    maxind2(2) = rec_cov_counter_age0devs;  
    maxind3(2) = nyr_recdevs;
    maxind2(3) = recruitment_covariate_counter;
    maxind2(4) = mortality_covariate_counter;
    maxind2(5) = mortality_covariate_counter;  
    maxind3(5) = nyr_tobefit;
    maxind2(6) = recruitment_covariate_counter;
    maxind2(7) = mortality_covariate_counter;
    maxind2(8) = rec_cov_counter_age0devs;
    maxind2(9) = mortality_covariate_counter;
    maxind2(10) = nyr_tobefit-1;

    if(vhs_rec_cotv==2){
        maxind2(11) = nyr_tobefit-1; // For time-varying - multiple recovery probabilities are estimated
    }  
    nyr_vhs = maxind2(11);

    maxind2(12) = nyr_tobefit-1;
    if(ich_rec_cotv==2){ 
        maxind2(13) = nyr_tobefit-1; // For time-varying - multiple recovery probabilities are estimated
    }  
    nyr_ich = maxind2(13);

 END_CALCS         

    // Containers for indicator on covariate variables
    vector beta_mortality_ind(1,mortality_covariate_counter)
    vector beta_recruit_ind(1,recruitment_covariate_counter)

    // Containers for array (vector and matrix) parameters
    3darray pars_2_init(1,npar2,minind2,maxind2,minind3,maxind3);
    ivector par_order(1,totpar)
    ivector par_type(1,totpar)

 LOCAL_CALCS
    // Setup 3D array specifying initial values for vector & matrix parameters
    // Initial values for each parameters are stored in a 3D array. The dimensions are:
    // 1 = Each Parameter
    // 2 = 1st dimension of parameter array (the vector, or rows of matrix)
    // 3 = 2nd dimension of parameter array (cols of matrix)
    j.fill_seqadd(1,1);
    for(int i=1; i<=npar2; i++){
        int itemp=1;
        for(int k=minind2(i); k<=maxind2(i); k++){
            for(int m=minind3(i); m<=maxind3(i); m++){
                pars_2_init(i,k,m) = pars_2(j(4),itemp);
                if(pars_2(j(4)).indexmax()>1){
                    itemp+=1;
                };
            }
        }
        // cout << "Filling 3D array " << i << endl << pars_2_init(i) << endl << endl;
        j+=4;
    }

    // Write PIN file if selected
    if(pin_write){
        // Vector of indices in pars_1 and pars_2_init to loop through for writing .PIN file
        // Order in PARAMETER_SECTION:
        //               1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45
        par_order.fill("{1,2,3,4,5,6,7,8,9,10,11,12,1 ,13,14,15,16,17,18,19,20,21,2 ,22,23,3 ,4 ,5 ,24,6 ,7 ,8 ,9 ,10,11,25,26,27,28,12,13,29,30,31,32}");
        par_type.fill( "{1,1,1,1,1,1,1,1,1,1 ,1 ,1 ,2 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,2 ,1 ,1 ,2 ,2 ,2 ,1 ,2 ,2 ,2 ,2 ,2 ,2 ,1 ,1 ,1 ,1 ,2 ,2 ,1 ,1 ,1 ,1 }");
        
        ofstream write_pin("PWS_ASA.PIN",ios::trunc);
        
        for(int i=1; i<=totpar; i++){
            int itemp = par_order(i);

            // If parameter is a number...
            if(par_type(i)==1){
                write_pin << pars_1(itemp,1) << endl;
            }else if(par_type(i)==2){
                // If parameter is a vector, write columnwise in .PIN
                if(pars_2_dim(itemp)==1){ 
                    write_pin << column(pars_2_init(itemp),1) << endl;
                // If parameter is a matrix, columns go columnwise, and rows rowwise
                }else if(pars_2_dim(itemp)==2){ 
                    dmatrix temp(minind2(itemp),maxind2(itemp),minind3(itemp),maxind3(itemp));
                    temp =  pars_2_init(itemp);
                    for(int k=minind2(itemp); k<=maxind2(itemp); k++){
                        write_pin << temp(k) << endl;
                    }
                }
                // cout << "Par " << itemp << " Dimensions: " << minind2(itemp) << " " << maxind2(itemp) << " " << minind3(itemp) << " " << maxind3(itemp) << endl;
            }
            // cout << "Par write complete index: " << i << endl;
        }
    }

    int j=1;
    for(int i=1; i<=n_age0_covs; i++){
        if(age0_turn_on(i)==1){
            beta_recruit_ind(j)=i;
            //beta_age0_PIN(j)=0.1;
            j+=1;
        }
    }

    j=1;
    for(int i=1; i<=n_mor_covs; i++){
        if(mor_turn_on(i)==1){
            beta_mortality_ind(j)=i;
            //beta_mortality_PIN(j)=0.1;
            j+=1;
        }
    }


    // Different maturity starting values and bounds by maturity_model_type 
    if(maturity_model_type==1 || maturity_model_type==2){
        // Direct proportions
        pars_1(5,1) = 0.60; pars_1(5,2) = 0.01; pars_1(5,3) = 0.90;
        pars_1(6,1) = 0.90; pars_1(6,2) = 0.30; pars_1(6,3) = 1.00;
        pars_1(7,1) = 0.60; pars_1(7,2) = 0.01; pars_1(7,3) = 0.90;
        pars_1(8,1) = 0.90; pars_1(8,2) = 0.30; pars_1(8,3) = 1.00;
    }else if(maturity_model_type==3){
        // For logistic function form of maturity function
        pars_1(5,1) = 3.00; pars_1(5,2) = 2.00; pars_1(5,3) = 5.00;
        pars_1(6,1) = 1.00; pars_1(6,2) = 0.05; pars_1(6,3) = 4.00;
        pars_1(7,1) = 0.00; pars_1(7,2) = 0.00; pars_1(7,3) = 1.00;
        pars_1(8,1) = 0.00; pars_1(8,2) = 0.00; pars_1(8,3) = 1.00;
    }else if(maturity_model_type==4){
        pars_1(5,1) = 0.6562465; pars_1(5,2) = 0.01; pars_1(5,3) = 0.90; pars_1(5,4) = -3;
        pars_1(6,1) = 0.8137432; pars_1(6,2) = 0.30; pars_1(6,3) = 1.00; pars_1(6,4) = -3;
        pars_1(7,1) = 0.5036630; pars_1(7,2) = 0.01; pars_1(7,3) = 0.90; pars_1(7,4) = -3;
        pars_1(8,1) = 0.3818182; pars_1(8,2) = 0.30; pars_1(8,3) = 1.00; pars_1(8,4) = -3;

        pars_1(11,4) = 4;
        pars_1(12,4) = 4;
    }

    if(DD_Mat==1){
        pars_1(5,2) =   0.00; pars_1(5,3) =  3.00;
        pars_1(6,2) = -12.00; pars_1(6,3) = -7.00;
        pars_1(7,2) =  -4.00; pars_1(7,3) = 4.00;
        pars_1(8,2) =   0.20; pars_1(8,3) = 1.00;
    }

 END_CALCS


PARAMETER_SECTION
    // To use bounds & phases from ctl file, need to individually assign them from
    // the design matrices (pars_1 & pars_2_ctl) to double type variables.
    // This needs to be done before assigning each parameters.
    // While repetitive, I don't know a way to write a function that returns multiple outputs separately.
    // Plus, users can directly see the assignments (which are identical...)


    // |---------------------------------------------------------------------------|
    // | NATURAL MORTALITY PARAMETERS
    // |---------------------------------------------------------------------------|
    !! int ii=1, jj=1; double lb, ub, ph; 

    //Age-related instantaneous mortality
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number VHSV_age3_4_mort_93(lb,ub,ph)

    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number ICH_age5_8_mort_93(lb,ub,ph)

    //Estimate baseline adult mortality in log-space - very non-informative prior
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number Z_0_8(lb,ub,ph)   
    
    //implements a constraint from ADF&G model: .25 <= S_9+ <= .95*S_5+, which 
    //than age3-8 natural mortality says that plus group mortality must be larger
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number Z_9(lb,ub,ph)
 
    // |---------------------------------------------------------------------------|
    // | MATURITY PARAMETERS
    // |---------------------------------------------------------------------------|
    //Maturity parameters of age 3 and 4 before (per1) and after 1997 (per2)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number mat_par_1(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number mat_par_2(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number mat_par_3(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number mat_par_4(lb,ub,ph)

    // |---------------------------------------------------------------------------|
    // | SELECTIVITY PARAMETERS
    // |---------------------------------------------------------------------------|
    //Seine Vulnerability parameters
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number seine_vuln_alpha(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number seine_vuln_beta(lb,ub,ph)

    //Survey Vulnerability parameters
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number survey_vuln_alpha(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number survey_vuln_beta(lb,ub,ph)

    //Initial abundance parameters (age-3 for all yrs, 1980 ages 4 and 5+)
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector loginit_pop(1,5,lb,ub,ph)

    // |---------------------------------------------------------------------------|
    // | SURVEY SCALAR & CV PARAMETERS
    // |---------------------------------------------------------------------------|
    // Egg deposition additional variance term
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number eggdep_add_var(lb,ub,ph) // In the likelihoods I use sqrt(), so bound must be
                                          // positive, since is non-differentiable at 0...

    //Milt coefficient
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number logmdm_c(lb,ub,ph)         // infinity if mdm_c goes to zero
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number milt_add_var(lb,ub,ph)

    //Hydroacoustic scalers and additional variance parameters
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number ADFG_hydro_q(lb,ub,ph) 
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number ADFG_hydro_add_var(lb,ub,ph)

    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number PWSSC_hydro_q(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number PWSSC_hydro_add_var(lb,ub,ph)

    // Aerial juvenile survey (incorporated 12/2019)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number log_juvenile_q(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number juvenile_overdispersion(lb,ub,ph)
 
    // |---------------------------------------------------------------------------|
    // | REC DEVIATES & COVARIATE EFFECTS ON MORTALITY
    // |---------------------------------------------------------------------------| 
    // Covariate effects on recruits (age 0) or mortality
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_matrix annual_age0devs(1,rec_cov_counter_age0devs,1,nyr_recdevs,lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number log_MeanAge0(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number sigma_age0devs(lb,ub,ph)
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector beta_age0(1,recruitment_covariate_counter,lb,ub,ph)

    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector beta_mortality(1,mortality_covariate_counter,lb,ub,ph)
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_matrix annual_mortdevs(1,mortality_covariate_counter,1,nyr_tobefit,lb,ub,ph) // Just estimate the deviates for 1992-1995
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number sigma_mortdevs(lb,ub,ph)

    // Basically sets 2nd time block for different effect estimate on covariates - need to improve....
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector beta_age0_offset(1,recruitment_covariate_counter,lb,ub,ph)
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector beta_mortality_offset(1,mortality_covariate_counter,lb,ub,ph)
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector sigma_age0covar(1,rec_cov_counter_age0devs,lb,ub,ph) // Weight of weighted Sum of Squares fit to recruitment indices
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector sigma_morcovar(1,mortality_covariate_counter,lb,ub,ph) // Weight of weighted SS fit to mortality indices

    // |---------------------------------------------------------------------------|
    // | AGE-STRUCTURED DISEASE EFFECTS BASED ON VHSV SEROPREVALENCE & ICH. HOF. INFECTION
    // |---------------------------------------------------------------------------|
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector vhsv_infection_prob(vhs_start_est,nyr_tobefit-1,lb,ub,ph)
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector vhsv_recovery_prob(vhs_start_est,nyr_vhs,lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number vhsv_infection_vuln_a50(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number vhsv_infection_vuln_a95(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number vhsv_sample_vuln_a50(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number vhsv_sample_vuln_a95(lb,ub,ph)

    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector ich_infection_prob(ich_start_est,nyr_tobefit-1,lb,ub,ph)
    !! lb = pars_2_ctl(jj,2); ub = pars_2_ctl(jj,3); ph = pars_2_ctl(jj,4); jj++;
    init_bounded_vector ich_recovery_prob(ich_start_est,nyr_ich,lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number ich_infection_vuln_a50(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++; 
    init_bounded_number ich_infection_vuln_a95(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number ich_sample_vuln_a50(lb,ub,ph)
    !! lb = pars_1(ii,2); ub = pars_1(ii,3); ph = pars_1(ii,4); ii++;
    init_bounded_number ich_sample_vuln_a95(lb,ub,ph)


    // |---------------------------------------------------------------------------|
    // | NUMERIC VARIABLES
    // |---------------------------------------------------------------------------|
    number S_0_2 			// Survival rate of ages 0-2
    number S_3_8			// Survival rate of ages 3-8
    number S_9			    // Survival rate of 9+ age group

    number mdm_c 
    number MDMtemp_2
    number HtempADFG_num

    number not_below_this

    //To count instances of incurring the posfun penalty
    number penCount

    number seine_llk        // purse-seine survey likelihood
    number spawner_llk      // spawning survey likelihood
    number mdm_llk          // mdm survey likelihood
    number eggdep_llk       // egg deposition survey likelihood
    number ADFG_hydro_llk   // ADFG hydroacoustic survey likelihood
    number PWSSC_hydro_llk  // PWSSC hydroacoustic survey likelihood
    number juvenile_llk     // Aerial juvenile survey likelihood (incorporated 12/2019)
    number vhsv_llk         // Seroprevalence likelihood
    number ich_llk          // Ichthyphonus infection likelihood

    number age0_devs_penllk
    number mort_devs_penllk
    number age0_covar_prior
    number mort_covar_prior

    //number mean_log_recruitment
    //number temp2MeanLgRec
    //number mean_log_recruitment
    number mean_log_rec
    number temp2MeanLgRec
    number meanLgRec
    number projected_prefishery_biomass //is an sdreport variable


    // Scalar for seroprevalence - scales down total immunity to account for
    // likely underestimation of immunity by current serological assays
    number q_immunity

    number priors

    // |---------------------------------------------------------------------------|
    // | CALCULATED VECTORS
    // |---------------------------------------------------------------------------|
    vector init_age_0(1,nyr_tobefit)       //Estimates of Age 0's (millions)
    vector forecast_winter_survival_effect(1,nage)
    vector forecast_survival_winter(1,nage)
    vector age0_effect(1,nyr_tobefit)

    vector seine_vuln(1,nage)                           //Vulnerability of fish to seine commercial fishery
    vector survey_vuln(1,nage)                          //Vulnerability of fish to seine gear used for surveys
    vector seine_avilable_pop(1,nyr_tobefit)            //Total Available population for Seine catch for each year
    vector seine_catch(1,nyr_tobefit)                   //Seine-Catch, Estimated total catch
    vector prefishery_biomass(1,nyr_tobefit)           //Pre-Fishery Run Biomass, mt
    vector prefishery_biomass_TONS(1,nyr_tobefit)      //Management reference point - is prefishery_biomass  (pre-fishery run biomass) in tons
    vector postfishery_spawn_biomass(1,nyr_tobefit)     //Total Post-Fishery Spawning Biomass
    vector prefishery_spawn_pop(1,nyr_tobefit)          //Total Pre-fishery spawning pop (mil)
    vector seine_biomass(1,nyr_tobefit)                 //Sum(Weight-at-age*Seine age comp) is like Seine biomass per year
    
    // TODO: What is this and how does it differ from pre_fish_biomass (mt)
    vector Early_biomass(1,nyr_tobefit)                 //Pre-Fishery Biomass by year

    // Aerial juvenile survey (incorporated 12/2019)
    vector juvenile_pred(1,nyr_tobefit)

    // These are model estimates for each of these survey values for each year
    // Difference between these estimated values and observed values constitutes
    // likelihood calculations. 
    vector MDM_est(1,nyr_tobefit)                   //Mile-days of milt
    vector EGG_est(1,nyr_tobefit)                   //Egg Deposition - rowsums of eggdep_age_comp
    vector EggAC_sum(1,nyr_tobefit)                 //totol egg dep (both male and female) by year
    vector ADFG_HYDRO_est(1,nyr_tobefit)              //ADFG Hydroacoustic Estimates
    vector PWSSC_HYDRO_est(1,nyr_tobefit)             //PWSSC Hydroacoustics Estimates

    vector mdm_residuals(1,nyr_tobefit)             //Mile-days of milt
    vector eggdep_residuals(1,nyr_tobefit)          //Egg deposition
    vector ADFG_hydro_residuals(1,nyr_tobefit)      // ADFG Hydroacoustics
    vector PWSSC_hydro_residuals(1,nyr_tobefit)     // PWSSC Hydroacoustics

    // Extra temporary variables for likelihood calculations.
    vector seine_temp2(1,nyr_tobefit)
    vector spawn_temp2(1,nyr_tobefit)
    vector seine_temp3(1,nyr_tobefit)
    vector spawn_temp3(1,nyr_tobefit)

    //Analytical Sigmas
    vector egg_sd(1,nyr_tobefit)
    vector PWSSC_sd(1,nyr_tobefit)

    //Likelihood components 
    vector mdm_llk_ind(1,nyr_tobefit)
    vector egg_llk_ind(1,nyr_tobefit)
    vector ADFG_hydro_llk_ind(1,nyr_tobefit)
    vector PWSSC_hydro_llk_ind(1,nyr_tobefit)
    vector juvenile_llk_ind(1,nyr_tobefit) // Aerial juvenile survey (incorporated 12/2019)
    
    // Variables and vectors for calculating projected final year biomass
    // vector tempWgt(1,nage)
    // vector avgWgt5Yr(1,nage)
    // vector projected_N_y_a(1,nage)
    // vector projected_Early_Sp_biomass(1,nage)

    vector Mean_Age0(1,nyr_tobefit) 
    vector init_pop(1,5)
    vector agggregate_annual_age0devs(1,nyr_tobefit)
    vector temp_age0(1,recruitment_covariate_counter)

    // Age-specific vulnerability to VHSV transmission (i.e. proportion mixing with sympatric pop)
    vector Vul_sympat(1,nage)
    vector Vul_vhs_survey(1,nage)

    vector inf_inc_sp(1,nyr_tobefit)
    vector vhsprev_sp(1,nyr_tobefit)
    vector fatal_sp(1,nyr_tobefit)

    vector inf_inc_age3(1,nyr_tobefit)
    vector vhsprev_age3(1,nyr_tobefit)
    vector fatal_age3(1,nyr_tobefit)  

    // Age-specific vulnerability to Ichthyophonus (i.e. proportion mixing with sympatric pop)
    vector Vul_ich_sympat(1,nage)
    vector Vul_ich_survey(1,nage)

    vector inf_inc_sp_ich(1,nyr_tobefit)
    vector ichprev_sp(1,nyr_tobefit)
    vector fatal_sp_ich(1,nyr_tobefit)

    vector inf_inc_age3_ich(1,nyr_tobefit)
    vector ichprev_age3(1,nyr_tobefit)
    vector fatal_age3_ich(1,nyr_tobefit)

    // Prior vector (stores individual values)
    vector prior_pars_1(1,npar)
    vector prior_pars_2(1,npar2)

    // |---------------------------------------------------------------------------|
    // | CALCULATED MATRICES
    // |---------------------------------------------------------------------------|
    //Mortality and survival rates
    matrix summer_mortality_effect(1,nyr_tobefit,1,nage)
    matrix winter_mortality_effect(1,nyr_tobefit,1,nage)
    matrix summer_survival(1,nyr_tobefit,1,nage)     //Half-year Survival between spring and fall
    matrix winter_survival(1,nyr_tobefit,1,nage)     //Half-year Survival between fall and spring

    //Age matrices
    matrix maturity(1,nyr_tobefit,1,nage)                           //Maturity 
    matrix N_y_a(1,nyr_tobefit,1,nage)                              //Pre-Fishery Abundance (millions)
    matrix vuln_pop_age(1,nyr_tobefit,1,nage)                       //Vulnerability_Seine*N_y_a, avail pop at age
    matrix seine_age_comp(1,nyr_tobefit,1,nage)                     //Seine Age Composition
    matrix N_spawners_age(1,nyr_tobefit,1,nage)                     //Spawning Population (millions)
    matrix prefishery_spawning_biomass_age(1,nyr_tobefit,1,nage)
    matrix spawning_biomass_age(1,nyr_tobefit,1,nage)               //Spawning Biomass at age
    matrix prefishery_spawning_pop_age(1,nyr_tobefit,1,nage)        //Pre-Fishery spawning population (millions)
    matrix spawning_age_comp(1,nyr_tobefit,1,nage)                  //Spawning Age Composition
    matrix seine_biomass_age(1,nyr_tobefit,1,nage)                  //Weight-at-age*Seine age comp, is like Seine biomass-at-age
    matrix eggdep_age_comp(1,nyr_tobefit,1,nage)                    //Egg dep by age, like Egg deposition Age-Comp  
    matrix prefishery_biomass_age(1,nyr_tobefit,1,nage)             //N_y_a*weight_at_age, Pre-Fishery biomass-at-age

    //Temporary matrices/vectors for calculating likelihood components
    //matrix SeACR(1,nyr_tobefit,1,nage)      //Seine age comp
    //matrix SpACR(1,nyr_tobefit,1,nage)      //Spawning age comp
    matrix seine_comp_residuals(1,nyr_tobefit,1,nage)
    matrix spawner_comp_residuals(1,nyr_tobefit,1,nage)

    matrix annual_mortdevs_byage(1,nyr_tobefit,1,nage)

    matrix maturity_unobs(1,nyr_tobefit,1,nage)      //Maturity of unobserved population

    // Age-specific matrices related to epidemiological stages
    matrix vhsv_susceptibility_age(1,nyr_tobefit,1,nage)
    matrix vhsv_immunity_age(1,nyr_tobefit,1,nage)
    matrix vhsv_survival_age(1,nyr_tobefit,1,nage)
    matrix vhsv_pred(1,nyr_tobefit,1,n_age2)

    matrix ich_susceptibility_age(1,nyr_tobefit,1,nage)
    matrix ich_chronic_infection_age(1,nyr_tobefit,1,nage)
    matrix ich_survival_age(1,nyr_tobefit,1,nage)
    matrix ich_pred(1,nyr_tobefit,1,n_age2)

    // |---------------------------------------------------------------------------|
    // | OBJECTIVE FUNCTION
    // |---------------------------------------------------------------------------|
    objective_function_value full_llk

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

    if(no_estimation){  
        Mean_Age0 = exp(log_MeanAge0);
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
        calc_priors();
        calc_nll_components();  

        struct stat buf;
        if(stat("rep_out",&buf)!=0){
            system("mkdir rep_out");
        }

        ofstream deterministic_run("rep_out/deterministic_run.rep",ios::trunc);

        // Aerial juvenile survey (incorporated 12/2019)
        deterministic_run << "# Posterior Probability" << endl;
        deterministic_run << seine_llk +spawner_llk +eggdep_llk +ADFG_hydro_llk +PWSSC_hydro_llk +mdm_llk +age0_devs_penllk +mort_devs_penllk +priors +juvenile_llk +vhsv_llk +ich_llk << endl << endl;

        // Aerial juvenile survey (incorporated 12/2019)
        deterministic_run << "# Likelihood Components & Priors" << endl;
        deterministic_run << "# seine_llk  spawner_llk  eggdep_llk  ADFG_hydro_llk  PWSSC_hydro_llk  mdm_llk  age0_devs_penllk  mort_devs_penllk  juvenile_llk  vhsv_llk  ich_llk" << endl;
        deterministic_run << seine_llk << "  " << spawner_llk << "  " << eggdep_llk << "  " << ADFG_hydro_llk << "  " << PWSSC_hydro_llk << "  " << mdm_llk << "  " << age0_devs_penllk << "  " << mort_devs_penllk << "  " << juvenile_llk <<  " " <<  vhsv_llk << " " <<  ich_llk << endl << endl;

        deterministic_run << "# Pre-fishery Spawning Biomass (metric tons)" << endl;
        deterministic_run << prefishery_biomass  << endl << endl;
        deterministic_run << "# Post-fishery Spawning Biomass (metric tons)" << endl;
        deterministic_run << postfishery_spawn_biomass << endl << endl;
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

    Mean_Age0 = exp(log_MeanAge0);
    mdm_c = exp(logmdm_c);
    init_pop = exp(loginit_pop);

    full_llk=0;
    penCount = 0;

    calc_naturalmortality();

    if(DD_Mat==0){
        calc_maturity();
    }

    calc_selectivity();
    calc_statevariables();

    calc_surveyvalues();
    calc_nll_components();
    calc_priors();

    // Objective function ADMB will minimize
    // Aerial juvenile survey (incorporated 12/2019)
    full_llk = seine_llk +
               spawner_llk +
               eggdep_llk +
               ADFG_hydro_llk +
               PWSSC_hydro_llk +
               mdm_llk +
               age0_devs_penllk +
               mort_devs_penllk +
               age0_covar_prior +
               mort_covar_prior +
               juvenile_llk +
               vhsv_llk +
               ich_llk +
               priors;

    SSB_final_year = prefishery_biomass(nyr_tobefit); // is projected year's Pre-fishery Run Biomass in metric tons

    // "ifMCEvalPhase" goes inside the PROCEDURE_SECTION,
    if(mceval_phase()){
        project_biomass(); // Project current year biomass using last year's data
        write_chain_results();
    }

    if(clean != 0){
        remove("fmin.log");
        remove("PWS_ASA.b01");
        remove("PWS_ASA.b02");
        remove("PWS_ASA.b03");
        remove("PWS_ASA.b04");
        remove("PWS_ASA.bar");
        remove("PWS_ASA.eva");
        remove("PWS_ASA.htp");
        remove("PWS_ASA.p01");
        remove("PWS_ASA.p02");
        remove("PWS_ASA.p03");
        remove("PWS_ASA.p04");
        remove("PWS_ASA.r01");
        remove("PWS_ASA.r02");
        remove("PWS_ASA.r03");
        remove("PWS_ASA.r04");
    }


FUNCTION void calc_naturalmortality()
    summer_mortality_effect.initialize();
    winter_mortality_effect.initialize();
    summer_survival.initialize();
    winter_survival.initialize();
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
    // ----------
    //Half-Year Survival (Matches Excel version - uses desease data and only estimates plus group mortality)
    //S_3_8=exp(-0.5*Z_0_8); //Z_0_8 is a read-in param 
  
    for (int i=1;i<=nyr_tobefit;i++){
        // Z_annual(i)=(M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(k)*mor_covariates(i,k)*covariate_effect_byage(j,k);
        for(int j=1;j<=nage;j++){
            for(int k=1;k<=mortality_covariate_counter;k++){
                if(beta_mortality_ind(k)==0){
                }else if(mor_covariates(i,beta_mortality_ind(k))==-9){
                    summer_mortality_effect(i,j) += 0;
                    winter_mortality_effect(i,j) += 0;
                }else if(mor_season(beta_mortality_ind(k))==1){
                    summer_mortality_effect(i,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
                }else if(mor_season(beta_mortality_ind(k))==2){
                    winter_mortality_effect(i,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
                }else if(mor_season(beta_mortality_ind(k))==3){
                    // Additional conditional because each model year starts with summer, ends with winter even though 
                    // data are input for the calendar year
                    if(i==nyr_tobefit){
                        summer_mortality_effect(i,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
                    }else{
                        summer_mortality_effect(i,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
                        winter_mortality_effect(i+1,j) += (M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*mor_covariates(i,beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
                    } 
                }
            }
        }
    }

    // Needed because calendar year for winter effect is in the second year-season of the model. In other words, user can input 
    // covariate effect info for the year being forecasted.

    for(int j=1;j<=nage;j++){
        for(int k=1;k<=mortality_covariate_counter;k++){
            if(beta_mortality_ind(k)==0){
            }else if(nyr_tobefit_winter_covariate(k)==-9){
                forecast_winter_survival_effect(j) += 0;
            }else if(mor_season(beta_mortality_ind(k))==2){
                forecast_winter_survival_effect(j) += (M_change(nyr_tobefit)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*nyr_tobefit_winter_covariate(beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
            }else if(mor_season(beta_mortality_ind(k))==3){
                forecast_winter_survival_effect(j) += (M_change(nyr_tobefit)*beta_mortality_offset(k)+beta_mortality(k))*mor_turn_on(beta_mortality_ind(k))*nyr_tobefit_winter_covariate(beta_mortality_ind(k))*covariate_effect_byage(j,beta_mortality_ind(k));
            }
        }
    }
  
  
    if(mortality_covariate_model==1){
        // Previous form where covariates are incoporated as fixed variables - changed 07/05/2019
        for (int i=1;i<=nyr_tobefit;i++){
            for (int j=1;j<=(nage-1);j++){
                if(i==13){
                    if((j>=4) && (j<=5)){
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)+VHSV_age3_4_mort_93));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)));

                    }else if(j>=6){
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)+ICH_age5_8_mort_93));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)));
                    }else{
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)));
                    }
                }else if(i==14){
                    if((j>=4) && (j<=5)){
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)+VHSV_age3_4_mort_93));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)+VHSV_age3_4_mort_93));
                    }else if(j>=6){
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)+ICH_age5_8_mort_93));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)+ICH_age5_8_mort_93));
                    }else{
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)));
                    }
                }else if(i==15){
                    if((j>=4) && (j<=5)){
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)+VHSV_age3_4_mort_93));
                    }else if(j>=6){
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)+ICH_age5_8_mort_93));
                    }else{
                        summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)));
                        winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)));
                    }
                }else{
                    summer_survival(i,j)=exp(-(0.5*Z_0_8+summer_mortality_effect(i,j)));
                    winter_survival(i,j)=exp(-(0.5*Z_0_8+winter_mortality_effect(i,j)));
                }
                
                // Calculate survival for ages 3-8 for years 1981 and above  
                dvariable pen_Sur_1=0.0;
                dvariable high_survival_penalty_1=1-summer_survival(i,j);
                summer_survival(i,j)=1-posfun(high_survival_penalty_1, 0.01, pen_Sur_1);
                if(summer_survival(i,j)>=1){
                    penCount+=1;
                }
                full_llk+=1000*pen_Sur_1;

                dvariable pen_Sur_2=0.0;
                dvariable high_survival_penalty_2=1-winter_survival(i,j);
                winter_survival(i,j)=1-posfun(high_survival_penalty_2, 0.01, pen_Sur_2);
                if(winter_survival(i,j)>=1){
                    penCount+=1;
                }
                full_llk+=1000*pen_Sur_2;
            }

            if(i==1){
                summer_survival(i,nage)=exp(-(0.5*Z_9+summer_mortality_effect(i,nage)));
                winter_survival(i,nage)=exp(-(0.5*Z_9+winter_mortality_effect(i,nage)));

                dvariable pen_Sur_3=0.0;
                dvariable high_survival_penalty_3=1-summer_survival(i,nage);
                summer_survival(i,nage)=1-posfun(high_survival_penalty_3, 0.01, pen_Sur_3);
                if(summer_survival(i,nage)>=1){
                    penCount+=1;
                }
                full_llk+=1000*pen_Sur_3;

                dvariable pen_Sur_4=0.0;
                dvariable high_survival_penalty_4=1-winter_survival(i,nage);
                winter_survival(i,nage)=1-posfun(high_survival_penalty_4, 0.01, pen_Sur_4);
                if(winter_survival(i,nage)>=1){
                    penCount+=1;
                }
                full_llk+=1000*pen_Sur_4;
            }else {
                summer_survival(i,nage)=summer_survival(i,nage-1)*summer_survival(i-1,nage)/summer_survival(i-1,nage-1); // Plus age group
                winter_survival(i,nage)=winter_survival(i,nage-1)*winter_survival(i-1,nage)/winter_survival(i-1,nage-1); // Plus age group
            }
        }

    }else if(mortality_covariate_model==2){
        for (int i=1;i<=nyr_tobefit;i++){
            for (int j=1;j<=(nage-1);j++){
                // I must do this because I only want to estimate deviates on mortality for ages and years to which I am fitting data
                for (int k=1;k<=mortality_covariate_counter;k++){
                    if(beta_mortality_ind(k)==0){
                    }else if(mor_covariates(i,beta_mortality_ind(k))!=-9){
                        annual_mortdevs_byage(i,j)+=(M_change(i)*beta_mortality_offset(k)+beta_mortality(k))*covariate_effect_byage(j,beta_mortality_ind(k))*annual_mortdevs(k,i);
                    }
                }
                summer_survival(i,j)=exp(-(0.5*Z_0_8+annual_mortdevs_byage(i,j)));
                winter_survival(i,j)=exp(-(0.5*Z_0_8+annual_mortdevs_byage(i,j)));

                // Calculate survival for ages 3-8 for years 1981 and above  
                dvariable pen_Sur_1=0.0;
                dvariable high_survival_penalty_1=1-summer_survival(i,j);
                summer_survival(i,j)=1-posfun(high_survival_penalty_1, 0.01, pen_Sur_1);
                if(summer_survival(i,j)>=1){
                    penCount+=1;
                }
                full_llk+=1000*pen_Sur_1;

                dvariable pen_Sur_2=0.0;
                dvariable high_survival_penalty_2=1-winter_survival(i,j);
                winter_survival(i,j)=1-posfun(high_survival_penalty_2, 0.01, pen_Sur_2);
                if(winter_survival(i,j)>=1){
                    penCount+=1;
                }
                full_llk+=1000*pen_Sur_2;
            }

            if(i==1){

                summer_survival(i,nage)=exp(-(0.5*Z_9+annual_mortdevs_byage(i,nage)));
                winter_survival(i,nage)=exp(-(0.5*Z_9+annual_mortdevs_byage(i,nage)));

                dvariable pen_Sur_3=0.0;
                dvariable high_survival_penalty_3=1-summer_survival(i,nage);
                summer_survival(i,nage)=1-posfun(high_survival_penalty_3, 0.01, pen_Sur_3);
                if(summer_survival(i,nage)>=1){
                    penCount+=1;
                }
                full_llk+=1000*pen_Sur_3;

                dvariable pen_Sur_4=0.0;
                dvariable high_survival_penalty_4=1-winter_survival(i,nage);
                winter_survival(i,nage)=1-posfun(high_survival_penalty_4, 0.01, pen_Sur_4);
                if(winter_survival(i,nage)>=1){
                    penCount+=1;
                }
                full_llk+=1000*pen_Sur_4;
            }else {
                summer_survival(i,nage)=summer_survival(i,nage-1)*summer_survival(i-1,nage)/summer_survival(i-1,nage-1); // Plus age group
                winter_survival(i,nage)=winter_survival(i,nage-1)*winter_survival(i-1,nage)/winter_survival(i-1,nage-1); // Plus age group
            }
        }     
    }


FUNCTION void calc_maturity()
    maturity.initialize();
    maturity_unobs.initialize();
    
    if(maturity_model_type<=1){
        for (int i=1;i<=nyr_tobefit;i++){
            maturity(i)(1,3)=0; 
            
            maturity(i,4)=mat_par_1*mat_par_2;
            maturity(i,5)=mat_par_2;

            maturity(i,6)=1;
            maturity(i)(7,nage)=1;
        }
    }else if(maturity_model_type==2){
        // Maturity values before 1997 - parameterization in Muradian et al. 2017
        // Change to 14 if pre 1994, 17 pre 1997, nyr_tobefit if single Early_biomass
        int yr_mat_change=17;

        for (int i=1;i<=nyr_tobefit;i++){
            maturity(i)(1,3)=0; // Ages 0-2 fully immature

            // Maturity levels switch after yr_mat_change (1997) 
            if(i<=yr_mat_change){
                maturity(i,4)=mat_par_1*mat_par_2;
                maturity(i,5)=mat_par_2;
            }else {
                maturity(i,4)=mat_par_3*mat_par_4;
                maturity(i,5)=mat_par_4;
            }
            maturity(i,6)=1;
            maturity(i)(7,nage)=1;
        }
    }else if(maturity_model_type==3){
        for (int i=1;i<=nyr_tobefit;i++){
            maturity(i)(1,3)=0; 

            // Logistic maturity function
            // mat_par_1 = age @ 50% mature
            // mat_par_2 = age @ 95% mature
            
            maturity(i,4)=1/(1+exp(-log(19)*(3-mat_par_1)/(mat_par_2)));
            maturity(i,5)=1/(1+exp(-log(19)*(4-mat_par_1)/(mat_par_2)));
            maturity(i,6)=1;
            maturity(i)(7,nage)=1;
        }
    }else if(maturity_model_type==4){
        for (int i=1;i<=nyr_tobefit;i++){
            maturity(i)(1,3)=0; 
            maturity(i,4)=mat_par_1*mat_par_2;
            maturity(i,5)=mat_par_2;
            maturity(i)(6,nage)=1;

            maturity_unobs(i)(1,3)=0; 
            maturity_unobs(i,4)=mat_par_3*mat_par_4;
            maturity_unobs(i,5)=mat_par_4;
            maturity_unobs(i)(6,nage)=1;
        }
    }
  

FUNCTION void calc_selectivity()
  
    //Gear Selectivity
    seine_vuln.initialize();
    seine_vuln(1,3)=0; 
    for (int j=4;j<=nage;j++) {
        seine_vuln(j)=1/(1+exp(-1.0*seine_vuln_beta*(j-1-seine_vuln_alpha)));
    }

    //Survey Selectivity
    survey_vuln.initialize();
    survey_vuln(1,3)=0; 
    for (int j=4;j<=nage;j++) {
        survey_vuln(j)=1/(1+exp(-1.0*survey_vuln_beta*(j-1-survey_vuln_alpha)));
    }

FUNCTION void calc_statevariables()
    N_y_a.initialize();
    vuln_pop_age.initialize();
    seine_avilable_pop.initialize();
    seine_age_comp.initialize();
    seine_biomass_age.initialize();
    seine_biomass.initialize();
    seine_catch.initialize();
    N_spawners_age.initialize();
    prefishery_spawning_biomass_age.initialize();
    prefishery_biomass .initialize();
    prefishery_biomass_TONS .initialize();
    spawning_biomass_age.initialize();
    postfishery_spawn_biomass.initialize();
    prefishery_spawning_pop_age.initialize();
    prefishery_spawn_pop.initialize();
    spawning_age_comp.initialize();
    prefishery_biomass_age.initialize();
    Early_biomass.initialize();
    age0_effect.initialize();
    init_age_0.initialize();
    agggregate_annual_age0devs.initialize();

    not_below_this = 0.01;

    // Age-specific vulnerability to sympatric transmission
    Vul_sympat.initialize();
    Vul_ich_sympat.initialize();
    for (int a=1;a<=nage;a++){
        // Vul_sympat(a) = 1/(1+exp(-log(19)*((a-1)-vhsv_infection_vuln_a50)/(vhsv_infection_vuln_a95-vhsv_infection_vuln_a50)));
        Vul_ich_sympat(a) = 1.0;
        Vul_sympat(a) = 1.0;
    }

    // All new born fish are 100% susceptible - COME BACK
    // vhsv_susceptibility_age.colfill(1,1.0);
    vhsv_susceptibility_age.initialize();
    vhsv_immunity_age.initialize();
    ich_susceptibility_age.initialize();
    ich_chronic_infection_age.initialize();
    vhsv_susceptibility_age = 1.0;
    vhsv_survival_age = 1.0;
    ich_susceptibility_age = 1.0;
    ich_survival_age = 1.0;

    // Time-varying OR constant recovery probability
    dvar_vector rec_prob_vhs2(vhs_start_est,nyr_tobefit-1);
    dvar_vector rec_prob_ich2(ich_start_est,nyr_tobefit-1);
    if(vhs_rec_cotv==1){
        rec_prob_vhs2 = vhsv_recovery_prob(vhs_start_est);  // When assuming constant, only uses first parameter in vector
    }else{
        rec_prob_vhs2 = vhsv_recovery_prob;
    }

    if(ich_rec_cotv==1){
        rec_prob_ich2 = ich_recovery_prob(ich_start_est);
    }else{
        rec_prob_ich2 = ich_recovery_prob;
    }

    // INITIAL CONDITIONS
    // init_age_0 is the number of recruits in year i
    // Sur_age0_2 is the survival rate experienced by init_age_0(i) from age 0.
    // In other words, the mortality the larvae in year i had undergone. Referred to as recruit because previous model started at age 3
    for(int k=1;k<=recruitment_covariate_counter;k++){
      if(beta_recruit_ind(k)==0){continue;}

      if(age0_covariates(1,beta_recruit_ind(k))==-9){
        age0_effect(1) += 0;
      }else{
        age0_effect(1) += age0_turn_on(beta_recruit_ind(k))*(R_change(1)*beta_age0_offset(k)+beta_age0(k))*age0_covariates(1,beta_recruit_ind(k));
        
        // This is only used if recruitment_covariate_model==2
        agggregate_annual_age0devs(1) += (R_change(1)*beta_age0_offset(k)+beta_age0(k))*annual_age0devs(k,1);
      }
    }

    if(recruitment_covariate_model==1){
      // Form below for effects incorporated as fixed variables
      init_age_0(1) = Mean_Age0(1)*exp((age0_effect(1)+annual_age0devs(1,1)-0.5*square(sigma_age0devs))); 
    }else if(recruitment_covariate_model==2){
      init_age_0(1) = Mean_Age0(1)*exp(agggregate_annual_age0devs(1)-0.5*square(sigma_age0devs)); 
    }

    //Fills in row 1 of pre-fishery abundance
    N_y_a(1,1)=init_age_0(1);
    --N_y_a(1)(2,6)=init_pop;   // Fill in first year N_y_a subvector with estimated initial ages 1-5 
    N_y_a(1)(7,nage)=0;      

    if(DD_Mat==1){
        maturity.initialize();
        maturity(1)(1,3)=0;
        maturity(1,4)=1/(1+exp(mat_par_1+exp(mat_par_2)*sum(N_y_a(1)(1,nage))));
        maturity(1,5)=maturity(1,4)+(1-maturity(1,4))/(1+exp(mat_par_3));
        maturity(1)(6,nage)=1;
    }
    // Calculate pre-fishery spawning biomass - nest elem_prod because function only accepts 2 args
    if(maturity_model_type!=4){
        prefishery_spawning_biomass_age(1)(1,nage)=elem_prod(elem_prod(maturity(1)(1,nage),N_y_a(1)(1,nage)),weight_at_age(1)(1,nage));
    }else{
        for(int j=1;j<=nage;j++){
            prefishery_spawning_biomass_age(1,j)=(survey_vuln(j)*maturity(1,j)+(1-survey_vuln(j))*maturity_unobs(1,j))*N_y_a(1,j)*weight_at_age(1,j);
        }
    }

  
    prefishery_biomass (1)=sum(prefishery_spawning_biomass_age(1)(1,nage));

    prefishery_biomass_TONS [1]=prefishery_biomass [1]*2204.62/2000; // pre-fishery run biomass in TONS
    
    //Fills in row 1 of Catch Age-Composition from Purse-Seine
    vuln_pop_age(1)(1,nage)=elem_prod(seine_vuln,N_y_a(1)(1,nage));
    seine_avilable_pop(1)=sum(vuln_pop_age(1)(1,nage)); 

    seine_age_comp(1)(1,nage)=vuln_pop_age(1)(1,nage)/seine_avilable_pop(1);
    seine_biomass_age(1)(1,nage)=elem_prod(seine_age_comp(1)(1,nage),weight_at_age(1)(1,nage));
    seine_biomass=sum(seine_biomass_age(1));

    if(maturity_model_type!=4){
        prefishery_spawning_pop_age(1)(1,nage)=elem_prod(maturity(1)(1,nage),N_y_a(1)(1,nage)); //numerator of Spawning Age-Comp(spawning_age_comp)
    }else{
        prefishery_spawning_pop_age(1)(1,nage)=elem_prod(survey_vuln(1,nage),N_y_a(1)(1,nage)); //numerator of Spawning Age-Comp(spawning_age_comp)
    }
    prefishery_spawn_pop(1)=sum(prefishery_spawning_pop_age(1)(1,nage));
    spawning_age_comp(1)(1,nage)=prefishery_spawning_pop_age(1)(1,nage)/prefishery_spawn_pop(1);  //fills in Spawning Age-Comp

    //Row 1 of total Seine catch in millions 
    seine_catch(1)=seine_yield(1)/seine_biomass(1);

    int count = 7;
    double gillnet_catch_sum;
    double pound_catch_sum;
    double food_bait_catch_sum;
        
    // Now generate naturally spawning pop as correct structure for first four years
    N_spawners_age(1)(1,nage)=0;
    spawning_biomass_age(1)(1,nage)=0;

    postfishery_spawn_biomass(1)=0;
    dvar_vector total_catch = seine_age_comp(1)(4,6)*seine_catch(1)+
                              gillnet_catch(1)(4,6)+
                              pk*pound_catch(1)(4,6);

    N_spawners_age(1)(4,6)=elem_prod(maturity(1)(4,6),N_y_a(1)(4,6)-total_catch);

    for(int j=4;j<=6;j++){
        dvariable pen4=0.0;
        N_spawners_age(1,j)=posfun(N_spawners_age(1,j), .01, pen4);
        if(N_spawners_age(1,j) <= not_below_this){
            penCount+=1;
        }
        full_llk+=1000*pen4;
    } 
  
    //spawning_biomass_age(1)(4,6)=elem_prod(elem_prod(N_spawners_age(1)(4,6),weight_at_age(1)(4,6)),1/maturity(1)(4,6));
    spawning_biomass_age(1)(4,6)=elem_prod(N_spawners_age(1)(4,6),weight_at_age(1)(4,6));
    
    postfishery_spawn_biomass(1)=sum(spawning_biomass_age(1)(4,6)); // postfishery_spawn_biomass=rowsum(spawning_biomass_age); //Total naturally spawning biomass
    int m=7;

    // FILL IN THE REMAINING YEARS
    dvariable total_spring_catch;
    for(int i=2;i<=nyr_tobefit;i++){
        // init_age_0 is the number of recruits in year i

        if(i<=(nyr_recdevs)){
            temp_age0 = column(annual_age0devs,i);
        }else{
            temp_age0 = 0;
        }

        for(int k=1;k<=rec_cov_counter_age0devs;k++){
            if(beta_recruit_ind(k)==0){continue;}
            
            if(age0_covariates(i,beta_recruit_ind(k))==-9){
                age0_effect(i) += 0;
            }else{
                age0_effect(i) += age0_turn_on(beta_recruit_ind(k))*(R_change(i)*beta_age0_offset(k)+beta_age0(k))*age0_covariates(i,beta_recruit_ind(k));
                agggregate_annual_age0devs(i) += (R_change(i)*beta_age0_offset(k)+beta_age0(k))*temp_age0(k);
            }
        }

        if(recruitment_covariate_model==1){
            init_age_0(i) = Mean_Age0(i)*exp((age0_effect(i)+temp_age0(1)-0.5*square((sigma_age0devs))));
        }else if(recruitment_covariate_model==2){
            init_age_0(i) = Mean_Age0(i)*exp(agggregate_annual_age0devs(i)-0.5*square(sigma_age0devs)); 
        }

        if(i<=5){                     //Fill in years 2:5 as plus group advances from 6+ to 9+
            N_y_a(i,1)=init_age_0(i);
            for(int j=2;j<=count-1;j++){
                N_y_a(i,j)=((N_y_a(i-1,j-1)-(seine_age_comp(i-1,j-1)*seine_catch(i-1)+gillnet_catch(i-1,j-1)+pk*pound_catch(i-1,j-1)))*summer_survival(i-1,j-1)-food_bait_catch(i-1,j-1))*winter_survival(i,j-1);
            }
            
            //The last age-class each year (incrementally by year beginning with age 6 in 1981) is a plus group
            gillnet_catch_sum = gillnet_catch(i-1,count-1);
            pound_catch_sum = pound_catch(i-1,count-1);
            food_bait_catch_sum = food_bait_catch(i-1,count-1);
            for(int k=count;k<=nage;k++){
                gillnet_catch_sum += gillnet_catch(i-1,k);
                pound_catch_sum += pound_catch(i-1,k);  
                food_bait_catch_sum += food_bait_catch(i-1,k); 
            }
            int j=count;

            dvariable seine_catch_sum = seine_age_comp(i-1,j-1)*seine_catch(i-1);
            total_spring_catch = seine_catch_sum+gillnet_catch_sum+pk*pound_catch_sum;

            // Ignore disease survival here because no disease data
            N_y_a(i,j)=(((N_y_a(i-1,j-1)-(total_spring_catch))*summer_survival(i-1,j-1))-food_bait_catch_sum)*winter_survival(i,j-1); 

            for(int j=1;j<=count;j++){
                dvariable pen1=0.0;
                N_y_a(i,j)=posfun(N_y_a(i,j), 0.01, pen1);
                if(N_y_a(i,j)<=not_below_this){
                    penCount+=1;
                }       
                full_llk+=1000*pen1;
                vuln_pop_age(i,j)=N_y_a(i,j)*seine_vuln(j);
            }

            seine_avilable_pop=rowsum(vuln_pop_age); 
            for(int j=1;j<=count;j++){
                seine_age_comp(i,j)=vuln_pop_age(i,j)/seine_avilable_pop(i);
                seine_biomass_age(i,j)=seine_age_comp(i,j)*weight_at_age(i,j); 
            }
            seine_biomass=rowsum(seine_biomass_age);
                        
            seine_catch(i)=seine_yield(i)/seine_biomass(i);
            count+=1;

        }else if((i>5) && (i<=13)){
            //Fills in the rest of each of the above matrices (N_y_a, Seine Catch age-comp and Total Seine catch) in 1 loop... Thank you, Hulson!
            //Below fills from 1985 to 1992
            N_y_a(i,1)=init_age_0(i);
            for(int j=2;j<=nage-1;j++){
                // Still no disease data, so no disease survival component (11/15/2020)
                total_spring_catch = seine_age_comp(i-1,j-1)*seine_catch(i-1)+
                                               gillnet_catch(i-1,j-1)+
                                               pk*pound_catch(i-1,j-1);
                N_y_a(i,j) = ((N_y_a(i-1,j-1)-(total_spring_catch))*summer_survival(i-1,j-1)-food_bait_catch(i-1,j-1))*winter_survival(i,j-1);
            }
            // Plus group
            int j=nage;
            total_spring_catch = seine_age_comp(i-1,j-1)*seine_catch(i-1)+
                                            gillnet_catch(i-1,j-1)+
                                            pk*pound_catch(i-1,j-1);
            N_y_a(i,j)=((N_y_a(i-1,j-1)-(total_spring_catch))*summer_survival(i-1,j-1)-food_bait_catch(i-1,j-1))*winter_survival(i,j-1) + ((N_y_a(i-1,j)-(seine_age_comp(i-1,j)*seine_catch(i-1)+gillnet_catch(i-1,j)+pk*pound_catch(i-1,j)))*summer_survival(i-1,j)-food_bait_catch(i-1,j))*winter_survival(i,j);
            
            for(int j=1;j<=nage;j++){
                dvariable pen2=0.0;
                N_y_a(i,j)=posfun(N_y_a(i,j), 0.01, pen2);
                if(N_y_a(i,j)<=not_below_this){
                    penCount+=1;
                }
                full_llk+=1000*pen2;
                vuln_pop_age(i,j)=N_y_a(i,j)*seine_vuln(j);
            } 

            seine_avilable_pop=rowsum(vuln_pop_age); 
            for(int j=1;j<=nage;j++){
                seine_age_comp(i,j)=vuln_pop_age(i,j)/seine_avilable_pop(i);
                seine_biomass_age(i,j)=seine_age_comp(i,j)*weight_at_age(i,j);
            }
            seine_biomass=rowsum(seine_biomass_age);

            seine_catch(i)=seine_yield(i)/seine_biomass(i);

        }else if(i>13){
            //Below fills from 1993 to nyr_tobefit - matches ADF&G Excel matrix exactly
            N_y_a(i,1)=init_age_0(i);
            for(int j=2;j<=nage-1;j++){

                if((i-1)>=vhs_start_est){
                    vhsv_survival_age(i-1,j-1) = (1 - Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1)) + Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1);
                    vhsv_immunity_age(i,j) = ( vhsv_immunity_age(i-1,j-1)+Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1) ) / ( (1 - Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1))+Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1) );
                    vhsv_susceptibility_age(i,j) = 1 - vhsv_immunity_age(i,j);
                }

                if((i-1)>=ich_start_est){
                    ich_survival_age(i-1,j-1) = (1 - Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1)) + Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1)*rec_prob_ich2(i-1);
                    ich_chronic_infection_age(i,j) = ( ich_chronic_infection_age(i-1,j-1)+Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1)*rec_prob_ich2(i-1) ) / ( (1 - Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1))+Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1)*rec_prob_ich2(i-1) );
                    ich_susceptibility_age(i,j) = 1 - ich_chronic_infection_age(i,j);
                }
                total_spring_catch = seine_age_comp(i-1,j-1)*seine_catch(i-1)+
                                     gillnet_catch(i-1,j-1)+
                                     pk*pound_catch(i-1,j-1);
                N_y_a(i,j)=((N_y_a(i-1,j-1)-(total_spring_catch))*summer_survival(i-1,j-1)-food_bait_catch(i-1,j-1))*winter_survival(i,j-1) * vhsv_survival_age(i-1,j-1) * ich_survival_age(i-1,j-1);

            }
            
            // Plus group
            int j=nage;
            
            dvariable Nya_0_8 = ((N_y_a(i-1,j-1)-(seine_age_comp(i-1,j-1)*seine_catch(i-1)+gillnet_catch(i-1,j-1)+pk*pound_catch(i-1,j-1)))*summer_survival(i-1,j-1)-food_bait_catch(i-1,j-1))*winter_survival(i,j-1);
            dvariable Nya_9 = ((N_y_a(i-1,j)-(seine_age_comp(i-1,j)*seine_catch(i-1)+gillnet_catch(i-1,j)+pk*pound_catch(i-1,j)))*summer_survival(i-1,j)-food_bait_catch(i-1,j))*winter_survival(i,j);

            if((i-1)>=vhs_start_est){
                vhsv_survival_age(i-1,j-1) = (1 - Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1)) + Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1);
                vhsv_survival_age(i-1,j) = (1 - Vul_sympat(j)*vhsv_susceptibility_age(i-1,j)*vhsv_infection_prob(i-1)) + Vul_sympat(j)*vhsv_susceptibility_age(i-1,j)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1);

                dvariable immune_vhs_prenage = ( vhsv_immunity_age(i-1,j-1)+Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1) ) / ( (1 - Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1))+Vul_sympat(j-1)*vhsv_susceptibility_age(i-1,j-1)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1) );
                dvariable immune_vhs_curnage = ( vhsv_immunity_age(i-1,j)+Vul_sympat(j)*vhsv_susceptibility_age(i-1,j)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1) ) / ( (1 - Vul_sympat(j)*vhsv_susceptibility_age(i-1,j)*vhsv_infection_prob(i-1))+Vul_sympat(j)*vhsv_susceptibility_age(i-1,j)*vhsv_infection_prob(i-1)*rec_prob_vhs2(i-1) );
                vhsv_immunity_age(i,j) = (immune_vhs_prenage * Nya_0_8 + immune_vhs_curnage * Nya_9)/(Nya_0_8 + Nya_9);
                vhsv_susceptibility_age(i,j) = 1 - vhsv_immunity_age(i,j);
            }

            if((i-1)>=ich_start_est){
                ich_survival_age(i-1,j-1) = (1 - Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1)) + Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1)*rec_prob_ich2(i-1);
                ich_survival_age(i-1,j) = (1 - Vul_ich_sympat(j)*ich_susceptibility_age(i-1,j)*ich_infection_prob(i-1)) + Vul_ich_sympat(j)*ich_susceptibility_age(i-1,j)*ich_infection_prob(i-1)*rec_prob_ich2(i-1);

                dvariable chrinf_ich_prenage = ( ich_chronic_infection_age(i-1,j-1)+Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1)*rec_prob_ich2(i-1) ) / ( (1 - Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1))+Vul_ich_sympat(j-1)*ich_susceptibility_age(i-1,j-1)*ich_infection_prob(i-1)*rec_prob_ich2(i-1) );
                dvariable chrinf_ich_curnage = ( ich_chronic_infection_age(i-1,j)+Vul_ich_sympat(j)*ich_susceptibility_age(i-1,j)*ich_infection_prob(i-1)*rec_prob_ich2(i-1) ) / ( (1 - Vul_ich_sympat(j)*ich_susceptibility_age(i-1,j)*ich_infection_prob(i-1))+Vul_ich_sympat(j)*ich_susceptibility_age(i-1,j)*ich_infection_prob(i-1)*rec_prob_ich2(i-1) );
                ich_chronic_infection_age(i,j) = (chrinf_ich_prenage * Nya_0_8 + chrinf_ich_curnage * Nya_9)/(Nya_0_8 + Nya_9);
                ich_susceptibility_age(i,j) = 1 - ich_chronic_infection_age(i,j);
            }

            N_y_a(i,j) = Nya_0_8 * vhsv_survival_age(i-1,j-1) * ich_survival_age(i-1,j-1) +  Nya_9 * vhsv_survival_age(i-1,j) * ich_survival_age(i-1,j);

            // Check for impossible values (negative numbers) & blow up likelihood if present
            for(int j=1;j<=nage;j++){
                dvariable pen3=0.0;
                N_y_a(i,j)=posfun(N_y_a(i,j), 0.01, pen3);
                if(N_y_a(i,j)<=not_below_this){
                    penCount+=1;
                }
                full_llk+=1000*pen3;
                vuln_pop_age(i,j)=N_y_a(i,j)*seine_vuln(j);
            }

            seine_avilable_pop=rowsum(vuln_pop_age); 
            for(int j=1;j<=nage;j++){
                seine_age_comp(i,j)=vuln_pop_age(i,j)/seine_avilable_pop(i);
                seine_biomass_age(i,j)=seine_age_comp(i,j)*weight_at_age(i,j);
            }
            seine_biomass=rowsum(seine_biomass_age);

            seine_catch(i)=seine_yield(i)/seine_biomass(i);
        }

        //Pre-Fishery Spawning Biomass or Pre-Fishery Run Biomass, mt
        if(DD_Mat==1){
            maturity(i)(1,3)=0;
            maturity(i,4)=1/(1+exp(mat_par_1+exp(mat_par_2)*sum(N_y_a(i)(1,nage))));
            maturity(i,5)=maturity(i,4)+(1-maturity(i,4))/(1+exp(mat_par_3));
            maturity(i)(6,nage)=1;
        }

        //Pre-Fishery Spawning Age-Composition
        prefishery_biomass (i) = 0;
        prefishery_spawn_pop(i) = 0;
        for(int j=1;j<=nage;j++){
            if(maturity_model_type!=4){
                prefishery_spawning_biomass_age(i,j)=maturity(i,j)*N_y_a(i,j)*weight_at_age(i,j);
                prefishery_spawning_pop_age(i,j)=maturity(i,j)*N_y_a(i,j); //numerator of Spawning Age-Comp(spawning_age_comp)
            }else{
                prefishery_spawning_biomass_age(i,j)=(survey_vuln(j)*maturity(i,j)+(1-survey_vuln(j))*maturity_unobs(i,j))*N_y_a(i,j)*weight_at_age(i,j);
                prefishery_spawning_pop_age(i,j)=survey_vuln(j)*N_y_a(i,j); //numerator of Spawning Age-Comp(spawning_age_comp)
            }
            prefishery_biomass (i)+= prefishery_spawning_biomass_age(i,j); //Spawning Stock Biomass; sum over ages the pre-fishery spawning biomass by year
            prefishery_spawn_pop(i)+=prefishery_spawning_pop_age(i,j);
        }
    
        prefishery_biomass_TONS[i]=prefishery_biomass [i]*2204.62/2000; // pre-fishery run biomass in TONS
        for(int j=1;j<=nage;j++){
            spawning_age_comp(i,j)=prefishery_spawning_pop_age(i,j)/prefishery_spawn_pop(i);  //fills in Spawning Age-Comp
        }

        // Post-First half year Fisheries Spawning Population Estimates, called Naturally Spawning Pop in Excel model
        // N_spawners_age Spawning Population (in millions) - mature abundance with spring catch and impound removed
        // Set up the remaining age-classes each year (incrementally by year beginning with age 5 in 1980) are zero
        // First set first 4 rows all to zero
        if(i<=4){
            for(int j=1;j<=nage;j++){
                N_spawners_age(i,j)=0;
                spawning_biomass_age(i,j)=0;
            }
            //m++;
        }else if(i>=5){
            // Now fill in the remaining rows regularly
            m=nage;
        }

        postfishery_spawn_biomass(i)=0;
        for(int j=4;j<=m;j++){
            dvariable total_catch = seine_age_comp(i,j)*seine_catch(i)+gillnet_catch(i,j)+pound_catch(i,j);
            N_spawners_age(i,j)=maturity(i,j)*(N_y_a(i,j)-(total_catch));
            
            dvariable pen4=0.0;
            N_spawners_age(i,j)=posfun(N_spawners_age(i,j), .01, pen4);
            if(N_spawners_age(i,j)<=not_below_this){
                penCount+=1;
            }
            full_llk+=1000*pen4;

            spawning_biomass_age(i,j)=N_spawners_age(i,j)*weight_at_age(i,j);
            //spawning_biomass_age(i,j)=N_spawners_age(i,j)*weight_at_age(i,j)/maturity(i,j);
            
            postfishery_spawn_biomass(i)+=spawning_biomass_age(i,j); // postfishery_spawn_biomass=rowsum(spawning_biomass_age); //Total naturally spawning biomass
        } 
        m++;
        
    }

    // Infection incidence, fatality, and immunity of spawning population
    inf_inc_sp.initialize();
    fatal_sp.initialize();
    vhsprev_sp.initialize();
    inf_inc_age3.initialize();
    fatal_age3.initialize();
    vhsprev_age3.initialize();

    dvar_vector temp_inc(1,nage);
    dvar_vector temp_fat(1,nage);
    dvar_vector temp_imm(1,nage);

    for(int i=vhs_start_est;i<nyr_tobefit;i++){
        for(int j=1;j<=nage;j++){
            temp_inc(j) = maturity(i,j)*N_y_a(i,j)*Vul_sympat(j)*vhsv_susceptibility_age(i,j)*vhsv_infection_prob(i);
            temp_fat(j) = maturity(i,j)*N_y_a(i,j)*Vul_sympat(j)*vhsv_susceptibility_age(i,j)*vhsv_infection_prob(i)*(1-rec_prob_vhs2(i));
            temp_imm(j) = maturity(i,j)*N_y_a(i,j)*vhsv_immunity_age(i,j);
        }
        inf_inc_sp(i) = sum(temp_inc)/sum(N_spawners_age(i)(1,nage));
        fatal_sp(i) = sum(temp_fat)/sum(N_spawners_age(i)(1,nage));
        vhsprev_sp(i) = sum(temp_imm)/sum(N_spawners_age(i)(1,nage));

        inf_inc_age3(i) = (maturity(i,4)*N_y_a(i,4)*Vul_sympat(4)*vhsv_susceptibility_age(i,4)*vhsv_infection_prob(i))/N_spawners_age(i,4);
        fatal_age3(i) = (maturity(i,4)*N_y_a(i,4)*Vul_sympat(4)*vhsv_susceptibility_age(i,4)*vhsv_infection_prob(i)*(1-rec_prob_vhs2(i)))/N_spawners_age(i,4);
        vhsprev_age3(i) = (maturity(i,4)*N_y_a(i,4)*vhsv_immunity_age(i,4))/N_spawners_age(i,4);
    }


    inf_inc_sp_ich.initialize();
    fatal_sp_ich.initialize();
    ichprev_sp.initialize();
    inf_inc_age3_ich.initialize();
    fatal_age3_ich.initialize();
    ichprev_age3.initialize();

    for(int i=ich_start_est;i<nyr_tobefit;i++){
        for(int j=1;j<=nage;j++){
            temp_inc(j) = maturity(i,j)*N_y_a(i,j)*Vul_ich_sympat(j)*ich_susceptibility_age(i,j)*ich_infection_prob(i);
            temp_fat(j) = maturity(i,j)*N_y_a(i,j)*Vul_ich_sympat(j)*ich_susceptibility_age(i,j)*ich_infection_prob(i)*(1-rec_prob_ich2(i));
            temp_imm(j) = maturity(i,j)*N_y_a(i,j)*ich_chronic_infection_age(i,j);
        }
        inf_inc_sp_ich(i) = sum(temp_inc)/sum(N_spawners_age(i)(1,nage));
        fatal_sp_ich(i) = sum(temp_fat)/sum(N_spawners_age(i)(1,nage));
        ichprev_sp(i) = sum(temp_imm)/sum(N_spawners_age(i)(1,nage));

        inf_inc_age3_ich(i) = (maturity(i,4)*N_y_a(i,4)*Vul_ich_sympat(4)*ich_susceptibility_age(i,4)*ich_infection_prob(i))/N_spawners_age(i,4);
        fatal_age3_ich(i) = (maturity(i,4)*N_y_a(i,4)*Vul_ich_sympat(4)*ich_susceptibility_age(i,4)*ich_infection_prob(i)*(1-rec_prob_ich2(i)))/N_spawners_age(i,4);
        ichprev_age3(i) = (maturity(i,4)*N_y_a(i,4)*ich_chronic_infection_age(i,4))/N_spawners_age(i,4);
    }

FUNCTION void calc_surveyvalues()
    MDM_est.initialize(); 
    eggdep_age_comp.initialize();
    EggAC_sum.initialize();
    EGG_est.initialize();
    prefishery_biomass_age.initialize();
    Early_biomass.initialize();  
    ADFG_HYDRO_est.initialize();
    PWSSC_HYDRO_est.initialize();
    juvenile_pred.initialize();
    vhsv_pred.initialize();
    ich_pred.initialize();

    //Egg deposition - this data set patchily exists for 10 out of the nyr_tobefit years
    for(int i=1;i<=nyr_tobefit;i++){
        for(int j=1;j<=nage;j++){
            dvariable eac = 0;
            if(egg_dep(i)!=-9){
                eac=N_spawners_age(i,j)*fecundity(i,j);
            }
            eggdep_age_comp(i,j)=eac;
        }
    }

    // Compute egg deposition survey, MDM survey, and aerial juvenile survey (12/2019) data  
    EggAC_sum=rowsum(eggdep_age_comp);
    for(int i=1;i<=nyr_tobefit;i++){
        EGG_est(i)=0.000001*female_spawners(i)*EggAC_sum(i);
        MDM_est(i)=(1-female_spawners(i))*postfishery_spawn_biomass(i)/mdm_c;
        juvenile_pred(i)=N_y_a(i,2)*exp(log_juvenile_q); 
    }

    // ADFG & PWSSC Hydroacoustic Survey Biomass 
    for(int i=1;i<=nyr_tobefit;i++){
        for(int j=1;j<=nage;j++){
            if(maturity_model_type!=4){
                //prefishery_biomass_age(i,j)=N_y_a(i,j)*weight_at_age(i,j);
                prefishery_biomass_age(i,j)=maturity(i,j)*N_y_a(i,j)*weight_at_age(i,j);
            }else{
                prefishery_biomass_age(i,j)=(survey_vuln(j)*maturity(i,j)+(1-survey_vuln(j))*maturity_unobs(i,j))*N_y_a(i,j)*weight_at_age(i,j);
            }
        }
        Early_biomass=rowsum(prefishery_biomass_age); // includes mature and non-mature fish
        ADFG_HYDRO_est(i)=Early_biomass(i)*exp(ADFG_hydro_q);
        PWSSC_HYDRO_est(i)=Early_biomass(i)*exp(PWSSC_hydro_q);
    }

    // Seroprevalence survey values (11/2020)
    Vul_vhs_survey.initialize();
    Vul_vhs_survey(1, nage) = 1.0;
    // for (int a=1;a<=nage;a++){
    //     // Vul_vhs_survey(a) = 1/(1+exp(-log(19)*((a-1)-vhsv_sample_vuln_a50)/(vhsv_sample_vuln_a95-vhsv_sample_vuln_a50)));
    //     // Vul_vhs_survey(a) = 1/(1+exp(-log(19)*(3-mat_par_1)/(mat_par_2-mat_par_1)));
    //     Vul_vhs_survey(a) = 1.0;
    // }

    int k = 1;


    q_immunity = 1.0;

    for(int i=vhs_start_est;i<=nyr_tobefit;i++){
        for(int j=1;j<=nage;j++){
            vhsv_pred(i,k) = Vul_vhs_survey(j)*N_y_a(i,j)*vhsv_immunity_age(i,j)*q_immunity; // This assumes vhsprevalence samples target all mature individuals (i.e. maturity = survey vulnerability)
            k+=1;
            vhsv_pred(i,k) = Vul_vhs_survey(j)*N_y_a(i,j)*(1-vhsv_immunity_age(i,j)*q_immunity);
            k+=1;
        } 
        k=1;
        vhsv_pred(i)(1,n_age2) = vhsv_pred(i)(1,n_age2)/sum(vhsv_pred(i)(1,n_age2));
    }


    // Ichthyophonus infection values (11/2021)
    Vul_ich_survey.initialize();
    Vul_ich_survey(1, nage) = 1/(1+exp(-log(19)*(3-mat_par_1)/(mat_par_2-mat_par_1)));
    // for (int a=1;a<=nage;a++){
    //     Vul_ich_survey(a) = 1/(1+exp(-log(19)*(3-mat_par_1)/(mat_par_2-mat_par_1)));
    // }

    k = 1;

    for(int i=ich_start_est;i<=nyr_tobefit;i++){
        for(int j=1;j<=nage;j++){
            // Since ages 7 and older are grouped in Maya's data, group here too
            if(j<8){
                ich_pred(i,k) = Vul_ich_survey(j)*N_y_a(i,j)*ich_chronic_infection_age(i,j); // This assumes seroprevalence samples target all mature individuals (i.e. maturity = survey vulnerability)
                k+=1;
                ich_pred(i,k) = Vul_ich_survey(j)*N_y_a(i,j)*(1-ich_chronic_infection_age(i,j));
                k+=1;
            }else{
                ich_pred(i,k) += Vul_ich_survey(j)*N_y_a(i,j)*ich_chronic_infection_age(i,j);
                ich_pred(i,k+1) += Vul_ich_survey(j)*N_y_a(i,j)*(1-ich_chronic_infection_age(i,j));
            }
        } 
        k=1;
        ich_pred(i)(1,n_age2) = ich_pred(i)(1,n_age2)/sum(ich_pred(i)(1,n_age2));
    }


FUNCTION void calc_nll_components()
    age0_devs_penllk.initialize();
    mort_devs_penllk.initialize();
    age0_covar_prior.initialize();
    mort_covar_prior.initialize();
    seine_comp_residuals.initialize();
    spawner_comp_residuals.initialize();
    mdm_llk_ind.initialize();
    egg_llk_ind.initialize();
    ADFG_hydro_llk_ind.initialize();
    PWSSC_hydro_llk_ind.initialize();
    juvenile_llk_ind.initialize();

    seine_llk.initialize();
    spawner_llk.initialize();
    mdm_llk.initialize();
    eggdep_llk.initialize();
    ADFG_hydro_llk.initialize();
    PWSSC_hydro_llk.initialize();
    juvenile_llk.initialize();
    vhsv_llk.initialize();
    ich_llk.initialize();

    // Remove penalized lik for unconstrained rec devs - 12/22/2019
    if(sigma_age0devs==0){
        age0_devs_penllk = 0;
    }else{
        for(int i=1; i<=nyr_recdevs; i++){
            if(recruitment_covariate_model==1){
                age0_devs_penllk += log(sigma_age0devs)+0.5*square(annual_age0devs(1,i))/square(sigma_age0devs);
            }else if(recruitment_covariate_model==2){
                for (int k=1;k<=recruitment_covariate_counter;k++){
                    if(beta_recruit_ind(k)==0){continue;}
                    
                    if(age0_covariates(i,beta_recruit_ind(k))!=-9){
                        // age0_devs_penllk += log(sigma_age0devs)+0.5*square(annual_age0devs(i))/square(sigma_age0devs);
                    age0_covar_prior += log(sigma_age0covar(k))+0.5*square(annual_age0devs(k,i)-age0_covariates(i,beta_recruit_ind(k)))/square(sigma_age0covar(k));
                    }
                }
                age0_devs_penllk += log(sigma_age0devs)+0.5*square(colsum(annual_age0devs)(i))/square(sigma_age0devs);
            }
        }
    }

    if(mortality_covariate_model==2){
        for(int i=1; i<=nyr_tobefit; i++){
            for (int k=1;k<=mortality_covariate_counter;k++){
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


    //Seine Age Composition - this data set is very patchy
    for(int i=1;i<=nyr_tobefit;i++){
        for(int j=1;j<=nage;j++){
            if(purse_seine_age_comp(i,j)<=0){
                seine_comp_residuals(i,j)=0;
            }else if(seine_age_comp(i,j)==0){
                seine_comp_residuals(i,j)=0;
            }else{
                seine_comp_residuals(i,j)=purse_seine_age_comp(i,j)*(log(seine_age_comp(i,j))-log(purse_seine_age_comp(i,j)));
            }
        }
    }
    
    seine_temp2=rowsum(seine_comp_residuals);
    for(int i=1;i<=nyr_tobefit;i++){
        seine_temp3(i)=seine_ess(i)*seine_temp2(i);
    }
    seine_llk=-sum(seine_temp3);
  
    //Spawning Age Composition
    for(int i=1;i<=nyr_tobefit;i++){
        for(int j=1;j<=nage;j++){
            if(spawner_age_comp(i,j)<=0){
                spawner_comp_residuals(i,j)=0;
            }else if(spawning_age_comp(i,j)==0){
                spawner_comp_residuals(i,j)=0;
            }else{
                spawner_comp_residuals(i,j)=spawner_age_comp(i,j)*(log(spawning_age_comp(i,j))-log(spawner_age_comp(i,j)));
            }
        }
    }

    spawn_temp2=rowsum(spawner_comp_residuals);
    for(int i=1;i<=nyr_tobefit;i++){
        spawn_temp3(i)=spawner_ess(i)*spawn_temp2(i);
    }
    spawner_llk=-sum(spawn_temp3);

    //Mile-days of milt Likelihood component
    for(int i=1;i<=nyr_tobefit;i++) {
        mdm_residuals(i)=log(MDM_est(i))-log(mile_days_milt(i));
        mdm_llk_ind(i)=(mdm_residuals(i)*mdm_residuals(i))/(2*milt_add_var*milt_add_var)+log(milt_add_var);
    }
    MDMtemp_2=norm2(mdm_residuals);
    //M_VAR=MDMtemp_2/nyr_tobefit+(milt_add_var*milt_add_var);
    mdm_llk=nyr_tobefit*log(milt_add_var)+(.5*MDMtemp_2/(milt_add_var*milt_add_var));

    //Egg Deposition Likelihood component
    eggdep_llk=0;
    for(int i=1;i<=nyr_tobefit;i++){
        egg_sd(i)=0;
        eggdep_residuals(i)=0;
        egg_llk_ind(i)=0;
    }
    for(int i=5;i<=5;i++){
        egg_sd(i)=sqrt((cv_egg(i)*cv_egg(i))+(eggdep_add_var*eggdep_add_var));
        eggdep_residuals(i)=(log(EGG_est(i))-log(egg_dep(i)));
        egg_llk_ind(i)=log(egg_sd(i))+(.5*eggdep_residuals(i)*eggdep_residuals(i)/(egg_sd(i)*egg_sd(i)));
        eggdep_llk+=egg_llk_ind(i);
    }
    for(int i=9;i<=13;i++){
        egg_sd(i)=sqrt((cv_egg(i)*cv_egg(i))+(eggdep_add_var*eggdep_add_var));
        eggdep_residuals(i)=(log(EGG_est(i))-log(egg_dep(i)));
        egg_llk_ind(i)=log(egg_sd(i))+(.5*eggdep_residuals(i)*eggdep_residuals(i)/(egg_sd(i)*egg_sd(i)));
        eggdep_llk+=egg_llk_ind(i);
    }
    for(int i=15;i<=18;i++){
        egg_sd(i)=sqrt((cv_egg(i)*cv_egg(i))+(eggdep_add_var*eggdep_add_var));
        eggdep_residuals(i)=(log(EGG_est(i))-log(egg_dep(i)));
        egg_llk_ind(i)=log(egg_sd(i))+(.5*eggdep_residuals(i)*eggdep_residuals(i)/(egg_sd(i)*egg_sd(i)));
        eggdep_llk+=egg_llk_ind(i);
    }
    // egg_llk_ind = egg_llk_ind*10;

    //ADFG Hydroacoustic Survey Biomass Likelihood component
    for(int i=1;i<=hydADFG_start-1;i++){
        ADFG_hydro_residuals(i)=0;
    }

    int N_hydADFG=0;
    for(int i=hydADFG_start;i<=hydADFG_start+4;i++){ 
        // hyd~_start variable holds index of first year of
        // survey depending on data source
        // UPDATED 07/23/2015 to reflect 2 additional years of missing data
        ADFG_hydro_residuals(i)=log(ADFG_hydro(i))-log(ADFG_HYDRO_est(i));
        ADFG_hydro_llk_ind(i)=(ADFG_hydro_residuals(i)*ADFG_hydro_residuals(i))/(ADFG_hydro_add_var*ADFG_hydro_add_var*2)+log(ADFG_hydro_add_var);
        N_hydADFG+=1;
    }
    for(int i=hydADFG_start+4+1;i<=nyr_tobefit;i++){
        ADFG_hydro_residuals(i)=0;    //missing final 5 years of ADF&G data as of 07/2015
        ADFG_hydro_llk_ind(i)=0;
    }
    HtempADFG_num=norm2(ADFG_hydro_residuals);
    ADFG_hydro_llk=(N_hydADFG)*log(ADFG_hydro_add_var)+(0.5*HtempADFG_num/(ADFG_hydro_add_var*ADFG_hydro_add_var));
    // minus 5 to account for the missing final 5 years of ADF&G data

    //PWSSC Hydroacoustic Survey Biomass Likelihood component
    PWSSC_hydro_llk=0;

    int N_hydPWWSC=0;
    for(int i=1;i<=nyr_tobefit;i++){ // hyd_start variable holds index of first year of survey depending on data source
        if(PWSSC_hydro_cv(i)==-9) {
            PWSSC_sd(i)=0;
            PWSSC_hydro_residuals(i)=0;
            PWSSC_hydro_llk_ind(i)=0;
        }
        else{
            PWSSC_sd(i)=sqrt((PWSSC_hydro_cv(i)*PWSSC_hydro_cv(i))+(PWSSC_hydro_add_var*PWSSC_hydro_add_var));
            PWSSC_hydro_residuals(i)=(log(PWSSC_hydro(i))-log(PWSSC_HYDRO_est(i)));
            PWSSC_hydro_llk_ind(i)=log(PWSSC_sd(i))+(.5*PWSSC_hydro_residuals(i)*PWSSC_hydro_residuals(i)/(PWSSC_sd(i)*PWSSC_sd(i)));
            PWSSC_hydro_llk+=PWSSC_hydro_llk_ind(i);
            N_hydPWWSC+=1;
        }
    }

    // Aerial juvenile survey (incorporated 12/2019)
    juvenile_llk=0;
    for(int i=1;i<=nyr_tobefit;i++) {
        if(juvenile_survey_index(i)!=-9) {
            juvenile_llk_ind(i)=dnbinom(juvenile_survey_index(i),juvenile_pred(i),juvenile_overdispersion);
            juvenile_llk+=juvenile_llk_ind(i);
        }
    }
  
    // Seroprevalence obs likelihood - binomial (may need to come back to) (01/2021)
    dvar_matrix vhs_llk_ind(1,nyr_tobefit,1,n_age2);
    dvar_matrix ich_llk_ind(1,nyr_tobefit,1,n_age2);

    for(int i=vhs_start;i<=nyr_tobefit;i++){
        // vhs_obs_comp(i)(1,n_age2) = vhs_obs(i)(1,n_age2)/vhs_ss(i);

        dvar_vector vhs_obs_comp(1,n_age2);
        vhs_obs_comp = vhs_obs(i)(1,n_age2);

        dvar_vector vhs_pred_comp(1,n_age2);
        vhs_pred_comp = vhsv_pred(i)(1,n_age2);

        dvariable x_bin;
        dvariable n_bin;
        dvariable p_bin;

        for(int j=1;j<=nage;j++){
            if(vhs_obs_comp(j*2-1)<0 | (vhs_pred_comp(j*2-1)==0 && vhs_obs_comp(j*2-1)==0)){
                vhs_llk_ind(i,j)=0;
            }else if(vhs_pred_comp(j*2-1)==0 && vhs_obs_comp(j*2-1)>0){
                if(pars_2_ctl(10,4)>0){
                    vhs_llk_ind(i,j)=1000;
                }else{
                    vhs_llk_ind(i,j)=0; // Means infection parameter is fixed at 0, so should not contribute to likelihood
                }
            }else{
                // Binomial
                x_bin = vhs_obs_comp(j*2-1)*vhsv_antibody_ess(i);
                n_bin = (vhs_obs_comp(j*2-1)+vhs_obs_comp(j*2))*vhsv_antibody_ess(i);
                p_bin = vhs_pred_comp(j*2-1)/(vhs_pred_comp(j*2-1) + vhs_pred_comp(j*2));
                vhs_llk_ind(i,j) = dbinom(x_bin,n_bin,p_bin);
            }
        }
        // cout << vhs_llk_ind << endl;
    }

    for(int i=ich_start;i<=nyr_tobefit;i++){
        dvar_vector ich_obs_comp(1,n_age2);
        ich_obs_comp = ich_obs(i)(1,n_age2);

        dvar_vector ich_pred_comp(1,n_age2);
        ich_pred_comp = ich_pred(i)(1,n_age2);

        dvariable x_bin;
        dvariable n_bin;
        dvariable p_bin;

        for(int j=1;j<=nage;j++){
            if(ich_obs_comp(j*2-1)<=0 | ich_pred_comp(j*2-1)==0){
                ich_llk_ind(i,j)=0;
            }else{
                // Binomial
                x_bin = ich_obs_comp(j*2-1)*ich_ess(i);
                n_bin = (ich_obs_comp(j*2-1)+ich_obs_comp(j*2))*ich_ess(i);
                p_bin = ich_pred_comp(j*2-1)/(ich_pred_comp(j*2-1) + ich_pred_comp(j*2));
                ich_llk_ind(i,j) = dbinom(x_bin,n_bin,p_bin);
            }
        }
    }

    // Binomial
    vhsv_llk = sum(rowsum(vhs_llk_ind));
    ich_llk = sum(rowsum(ich_llk_ind));

FUNCTION void calc_priors()
    priors.initialize();
    prior_pars_1.initialize();
    prior_pars_2.initialize();

    // Single valued parameters
    prior_pars_1(1) = prior_dist(pars_1(1 ,5), VHSV_age3_4_mort_93, pars_1(1 ,6), pars_1(1 ,7));
    prior_pars_1(2) = prior_dist(pars_1(2 ,5), ICH_age5_8_mort_93,  pars_1(2 ,6), pars_1(2 ,7));
    prior_pars_1(3) = prior_dist(pars_1(3 ,5), Z_0_8,               pars_1(3 ,6), pars_1(3 ,7));
    prior_pars_1(4) = prior_dist(pars_1(4 ,5), Z_9,                 pars_1(4 ,6), pars_1(4 ,7));
    prior_pars_1(5) = prior_dist(pars_1(5 ,5), mat_par_1,           pars_1(5 ,6), pars_1(5 ,7));
    prior_pars_1(6) = prior_dist(pars_1(6 ,5), mat_par_2,           pars_1(6 ,6), pars_1(6 ,7));

    if(maturity_model_type==2){
        prior_pars_1(7) = prior_dist(pars_1(7,5), mat_par_3, pars_1(7,6), pars_1(7,7));
        prior_pars_1(8) = prior_dist(pars_1(8,5), mat_par_4, pars_1(8,6), pars_1(8,7));
    }

    prior_pars_1(9)  = prior_dist(pars_1(9 ,5), seine_vuln_alpha,           pars_1(9 ,6), pars_1(9 ,7));
    prior_pars_1(10) = prior_dist(pars_1(10,5), seine_vuln_beta,            pars_1(10,6), pars_1(10,7));
    prior_pars_1(11) = prior_dist(pars_1(11,5), survey_vuln_alpha,          pars_1(11,6), pars_1(11,7));
    prior_pars_1(12) = prior_dist(pars_1(12,5), survey_vuln_beta,           pars_1(12,6), pars_1(12,7));
    prior_pars_1(13) = prior_dist(pars_1(13,5), eggdep_add_var,             pars_1(13,6), pars_1(13,7));
    prior_pars_1(14) = prior_dist(pars_1(14,5), logmdm_c,                   pars_1(14,6), pars_1(14,7));
    prior_pars_1(15) = prior_dist(pars_1(15,5), milt_add_var,               pars_1(15,6), pars_1(15,7));
    prior_pars_1(16) = prior_dist(pars_1(16,5), ADFG_hydro_q,               pars_1(16,6), pars_1(16,7));
    prior_pars_1(17) = prior_dist(pars_1(17,5), ADFG_hydro_add_var,         pars_1(17,6), pars_1(17,7));
    prior_pars_1(18) = prior_dist(pars_1(18,5), PWSSC_hydro_q,              pars_1(18,6), pars_1(18,7));
    prior_pars_1(19) = prior_dist(pars_1(19,5), PWSSC_hydro_add_var,        pars_1(19,6), pars_1(19,7));
    prior_pars_1(20) = prior_dist(pars_1(20,5), log_juvenile_q,             pars_1(20,6), pars_1(20,7));
    prior_pars_1(21) = prior_dist(pars_1(21,5), juvenile_overdispersion,    pars_1(21,6), pars_1(21,7));
    prior_pars_1(22) = prior_dist(pars_1(22,5), log_MeanAge0,               pars_1(22,6), pars_1(22,7));
    prior_pars_1(23) = prior_dist(pars_1(23,5), sigma_age0devs,             pars_1(23,6), pars_1(23,7));
    prior_pars_1(24) = prior_dist(pars_1(24,5), sigma_mortdevs,             pars_1(24,6), pars_1(24,7));
    prior_pars_1(25) = prior_dist(pars_1(25,5), vhsv_infection_vuln_a50,    pars_1(25,6), pars_1(25,7));
    prior_pars_1(26) = prior_dist(pars_1(26,5), vhsv_infection_vuln_a95,    pars_1(26,6), pars_1(26,7));
    prior_pars_1(27) = prior_dist(pars_1(27,5), vhsv_sample_vuln_a50,       pars_1(27,6), pars_1(27,7));
    prior_pars_1(28) = prior_dist(pars_1(28,5), vhsv_sample_vuln_a95,       pars_1(28,6), pars_1(28,7));
    prior_pars_1(29) = prior_dist(pars_1(29,5), ich_infection_vuln_a50,     pars_1(29,6), pars_1(29,7));
    prior_pars_1(30) = prior_dist(pars_1(30,5), ich_infection_vuln_a95,     pars_1(30,6), pars_1(30,7));
    prior_pars_1(31) = prior_dist(pars_1(31,5), ich_sample_vuln_a50,        pars_1(31,6), pars_1(31,7));
    prior_pars_1(32) = prior_dist(pars_1(32,5), ich_sample_vuln_a95,        pars_1(32,6), pars_1(32,7));
    
    dvariable partemp;
    partemp.initialize();

    int par_i;

    // loginit_pop
    par_i=1;
    for(int i=1;i<=5;i++){
        partemp = loginit_pop(i);
        prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
    }

    // annual_age0devs
    par_i=2;
    for(int i=1;i<=rec_cov_counter_age0devs;i++){
        for(int j=1;j<=(nyr_tobefit-3);j++){
            partemp = annual_age0devs(i,j);
            prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
        }
    }

    // beta_age0
    par_i=3;
    for(int i=1;i<=recruitment_covariate_counter;i++){
        partemp = beta_age0(i);
        prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
    }

    // beta_mortality
    par_i=4;
    for(int i=1;i<=mortality_covariate_counter;i++){
        partemp = beta_mortality(i);
        prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
    }

    // annual_mortdevs
    par_i=5;
    for(int i=1;i<=mortality_covariate_counter;i++){
        for(int j=1;j<=(nyr_tobefit);j++){
            partemp = annual_mortdevs(i,j);
            prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
        }
    }

    // beta_age0_offset
    par_i=6;
    for(int i=1;i<=recruitment_covariate_counter;i++){
        partemp = beta_age0_offset(i);
        prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
    }

    // beta_mortality_offset
    par_i=7;
    for(int i=1;i<=mortality_covariate_counter;i++){
        partemp = beta_mortality_offset(i);
        prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
    }

    // sigma_age0covar
    par_i=8;
    for(int i=1;i<=rec_cov_counter_age0devs;i++){
        partemp = sigma_age0covar(i);
        prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
    }

    // sigma_morcovar
    par_i=9;
    for(int i=1;i<=mortality_covariate_counter;i++){
        partemp = sigma_morcovar(i);
        prior_pars_2(par_i) += prior_dist(pars_2_ctl(par_i,5), partemp, pars_2_ctl(par_i,6), pars_2_ctl(par_i,7));
    }

    // VHS infection probabilities (each year)
    for(int i=vhs_start_est;i<nyr_tobefit;i++){
        partemp = vhsv_infection_prob(i);
        prior_pars_2(10) += prior_dist(pars_2_ctl(10,5), partemp, pars_2_ctl(10,6), pars_2_ctl(10,7));
    }

    // VHS recovery probability
    if(vhs_rec_cotv==1){
        partemp = vhsv_recovery_prob(vhs_start_est);
        prior_pars_2(11) += prior_dist(pars_2_ctl(11,5), partemp, pars_2_ctl(11,6), pars_2_ctl(11,7));
    }else{
        for(int i=vhs_start_est;i<nyr_tobefit;i++){
            partemp = vhsv_recovery_prob(i);
            prior_pars_2(11) += prior_dist(pars_2_ctl(11,5), partemp, pars_2_ctl(11,6), pars_2_ctl(11,7));
        }
    }

    // I. hoferi infection probabilities (each year)
    for(int i=ich_start_est;i<nyr_tobefit;i++){
        partemp = ich_infection_prob(i);
        prior_pars_2(12) += prior_dist(pars_2_ctl(12,5), partemp, pars_2_ctl(12,6), pars_2_ctl(12,7));
    }

    // I. hoferi recovery probability
    if(ich_rec_cotv==1){
        partemp = ich_recovery_prob(ich_start_est);
        prior_pars_2(13) += prior_dist(pars_2_ctl(13,5), partemp, pars_2_ctl(13,6), pars_2_ctl(13,7));
    }else{
        for(int i=ich_start_est;i<nyr_tobefit;i++){
            partemp = ich_recovery_prob(i);
            prior_pars_2(13) += prior_dist(pars_2_ctl(13,5), partemp, pars_2_ctl(13,6), pars_2_ctl(13,7));
        }
    }
    
    priors = sum(prior_pars_1) + sum(prior_pars_2);

  
FUNCTION dvariable prior_dist(const int& priortype, const prevariable& priorval, const double& priorp1, const double& priorp2)
    dvariable ptmp;
    ptmp.initialize();
    // match up prior type designated in the control file with distribution
    switch(priortype){
        case 1:   //normal
            ptmp = dnorm(priorval,priorp1,priorp2);
            break;
            
        case 2:   //lognormal CHANGED RF found an error in dlnorm prior. rev 116
            ptmp = dlnorm(priorval,priorp1,priorp2);
            break;
            
        case 3:   //beta distribution (assumed already to be on 0-1 scale)
            ptmp = dbeta(priorval,priorp1,priorp2);
            break;
            
        case 4:   //gamma distribution
            ptmp = dgamma(priorval,priorp1,priorp2);
            break;
    }
    return(ptmp);

FUNCTION project_biomass    
    //Use the average differences across the last 5 years...
    dvector tempWgt(1,nage);
    dvector avgWgt5Yr(1,nage);
    dvar_vector projected_N_y_a(1, nage);
    dvar_vector projected_Early_Sp_biomass(1, nage);

    for (int a=1; a<=nage; a++){
        tempWgt(a) = 0;
        for (int i=0; i<=4; i++){
            tempWgt(a) += weight_at_age(nyr_tobefit-i, a);
        }
        avgWgt5Yr(a) = tempWgt(a)/5;
    }
    
    // Now using mean of log-recruits for projection; this should be near the median if the recruits are log-Normally distributed.
    //cout << init_age_0(33)<< endl;
    dvariable mean_log_rec = 1;
    for(int i=nyr_tobefit-9;i<=nyr_tobefit;i++){
        mean_log_rec *= N_y_a(i,4);
    }
    mean_log_rec = log(mean_log_rec)/10;
    mean_log_rec = exp(temp2MeanLgRec);
    //cout << endl << "meanLgRec: " << meanLgRec << " "<< endl;

    // Plug age-3 biomass into first element of N_y_a vector for troubleshooting and display in report file
    projected_N_y_a(4) = meanLgRec;
    
    // Plug age-3 biomass into first element of vector for calculations
    projected_Early_Sp_biomass(4) = maturity(nyr_tobefit,4)*meanLgRec*avgWgt5Yr(4);
    
    // Fill in half-year survival rates for forecast year
    for(int j=1;j<=nage;j++){
        forecast_survival_winter(j)=exp(-(0.5*Z_0_8+forecast_winter_survival_effect(j)));
    }

    // Call numbers this year using last year's info for ages 4-8
    for(int j=5;j<=nage-1;j++){
        projected_N_y_a(j)=((N_y_a(nyr_tobefit,j-1)-(seine_age_comp(nyr_tobefit,j-1)*seine_catch(nyr_tobefit)+gillnet_catch(nyr_tobefit,j-1)+pk*pound_catch(nyr_tobefit,j-1)))*summer_survival(nyr_tobefit,j-1)-food_bait_catch(nyr_tobefit,j-1))*forecast_survival_winter(j-1);
    }
    
    // Calc numbers this year using last year's info for plus group
    for(int j=nage;j<=nage;j++){
        projected_N_y_a(j)=((N_y_a(nyr_tobefit,j-1)-(seine_age_comp(nyr_tobefit,j-1)*seine_catch(nyr_tobefit)+gillnet_catch(nyr_tobefit,j-1)+pk*pound_catch(nyr_tobefit,j-1)))*summer_survival(nyr_tobefit,j-1)-food_bait_catch(nyr_tobefit,j-1))*forecast_survival_winter(j-1)+((N_y_a(nyr_tobefit,j)-(seine_age_comp(nyr_tobefit,j)*seine_catch(nyr_tobefit)+gillnet_catch(nyr_tobefit,j)+pk*pound_catch(nyr_tobefit,j)))*summer_survival(nyr_tobefit,j)-food_bait_catch(nyr_tobefit,j))*forecast_survival_winter(j);
    }
    
    // Make it biomass
    for(int j=5;j<=nage;j++){
        projected_Early_Sp_biomass(j) = maturity(nyr_tobefit,j)*projected_N_y_a(j)*avgWgt5Yr(j);
    }
        
    // Take total pre-fishery biomass for projection year
    projected_prefishery_biomass = sum(projected_Early_Sp_biomass);


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
    ofstream MCMCreport13("mcmc_out/Ich_pred.csv",ios::app); 
    ofstream LLikReport("mcmc_out/llikcomponents.csv",ios::app);
    ofstream PriorReport("mcmc_out/priordensities.csv",ios::app);
    ofstream PFRReport("mcmc_out/PFRBiomass.csv",ios::app);
    ofstream indiv_LLikReport("mcmc_out/llik_observations.csv",ios::app);
    ofstream recruit_effect_report("mcmc_out/recruitment_effects.csv",ios::app);
    ofstream summer_survival_report("mcmc_out/adult_survival_effects_summer.csv",ios::app);
    ofstream winter_survival_report("mcmc_out/adult_survival_effects_winter.csv",ios::app);
    ofstream vhs_report1("mcmc_out/vhs_survival.csv",ios::app);
    ofstream vhs_report2("mcmc_out/vhs_infection.csv",ios::app);
    ofstream vhs_report3("mcmc_out/vhs_prev.csv",ios::app);
    ofstream vhs_report4("mcmc_out/vhs_fatality.csv",ios::app);
    ofstream vhs_report5("mcmc_out/vhs_infection_age3.csv",ios::app);
    ofstream vhs_report6("mcmc_out/vhs_prev_age3.csv",ios::app);
    ofstream vhs_report7("mcmc_out/vhs_fatality_age3.csv",ios::app);
    ofstream ich_report1("mcmc_out/ich_survival.csv",ios::app);
    ofstream ich_report2("mcmc_out/ich_infection.csv",ios::app);
    ofstream ich_report3("mcmc_out/ich_prev.csv",ios::app);
    ofstream ich_report4("mcmc_out/ich_fatality.csv",ios::app);
    ofstream ich_report5("mcmc_out/ich_infection_age3.csv",ios::app);
    ofstream ich_report6("mcmc_out/ich_prev_age3.csv",ios::app);
    ofstream ich_report7("mcmc_out/ich_fatality_age3.csv",ios::app);

    MCMCreport1 << milt_add_var << "," << eggdep_add_var << "," << ADFG_hydro_add_var << ","  <<  PWSSC_hydro_add_var << endl;
    
    // These files have values written by year across columns (and age in some cases)
    for (int i=1; i<nyr_tobefit; i++){
      MCMCreport2 << N_y_a(i,4) << ",";
      MCMCreport3 << ADFG_HYDRO_est(i) << ",";
      MCMCreport4 << PWSSC_HYDRO_est(i) << ",";
      MCMCreport5 << EGG_est(i) << ",";
      MCMCreport6 << MDM_est(i) << ",";
      MCMCreport7 << postfishery_spawn_biomass(i) << ",";

      for (int j=1; j<=nage; j++){
        MCMCreport8 << seine_age_comp(i,j) << ",";
        MCMCreport9 << spawning_age_comp(i,j) << ",";
        MCMCreport11 << N_y_a(i,j) << ",";
      }
      MCMCreport10 << juvenile_pred(i) << ",";

      for (int k=1; k<=n_age2; k++){
          MCMCreport12 << vhsv_pred(i,k) << ",";
          MCMCreport13 << ich_pred(i,k) << ",";
      }

      recruit_effect_report << exp(age0_effect(i))*Mean_Age0(i) << ",";
      for (int j=1; j<=nage; j++){
        summer_survival_report << summer_survival(i,j) << ",";
        winter_survival_report << winter_survival(i,j) << ",";
      }
    }

    // File in final year separately with 'endl' so next MCMC sample is written on the next line
    MCMCreport2 << N_y_a(nyr_tobefit,4) << endl; // this is the projected recruitment for the latest year
    MCMCreport3 << ADFG_HYDRO_est(nyr_tobefit) << endl;
    MCMCreport4 << PWSSC_HYDRO_est(nyr_tobefit) << endl;
    MCMCreport5 << EGG_est(nyr_tobefit) << endl;
    MCMCreport6 << MDM_est(nyr_tobefit) << endl;
    MCMCreport7 << postfishery_spawn_biomass(nyr_tobefit) << endl;

    for (int j=1; j<nage; j++){
        MCMCreport8 << seine_age_comp(nyr_tobefit,j) << ",";
        MCMCreport9 << spawning_age_comp(nyr_tobefit,j) << ",";
        MCMCreport11 << N_y_a(nyr_tobefit,j) << ","; 
    }

    MCMCreport8 << seine_age_comp(nyr_tobefit,nage) << endl;
    MCMCreport9 << spawning_age_comp(nyr_tobefit,nage) << endl;
    MCMCreport10 << juvenile_pred(nyr_tobefit) << endl;
    MCMCreport11 << N_y_a(nyr_tobefit,nage) << endl;

    for (int k=1; k<n_age2; k++){
        MCMCreport12 << vhsv_pred(nyr_tobefit,k) << ",";
        MCMCreport13 << ich_pred(nyr_tobefit,k) << ",";
    }
    
    MCMCreport12 << vhsv_pred(nyr_tobefit,n_age2) << endl;
    MCMCreport13 << ich_pred(nyr_tobefit,n_age2) << endl;

    recruit_effect_report << exp(age0_effect(nyr_tobefit))*Mean_Age0(nyr_tobefit) << endl;
    for (int j=1; j<=nage-1; j++){
      summer_survival_report << summer_survival(nyr_tobefit,j) << ",";
      winter_survival_report << winter_survival(nyr_tobefit,j) << ",";
    }
    summer_survival_report << summer_survival(nyr_tobefit,nage) << endl;
    winter_survival_report << winter_survival(nyr_tobefit,nage) << endl;


    // Now output the loglikelihood components
    LLikReport << -seine_llk << "," << -spawner_llk << "," << eggdep_llk << "," << ADFG_hydro_llk << "," << PWSSC_hydro_llk << "," << mdm_llk << ",";
    // LLikReport << age0_devs_penllk << "," << mort_devs_penllk << ",";
    // Aerial juvenile survey (12/2019) and seroprevalence samples (11/2020)
    LLikReport << juvenile_llk << "," << vhsv_llk << "," << ich_llk << "," << full_llk << endl;

    // Prior density reports
    for(int j=1; j<=npar; j++){
      PriorReport << prior_pars_1(j) << ",";
    }
    for(int j=1; j<=npar2; j++){
      PriorReport << prior_pars_2(j) << ",";
    }

    //indiv_LLikReport << age0_devs_penllk << "," << mort_devs_penllk << "," << Z_prior << "," << hydADFG_add_prior << "," << hydPWSSC_add_prior << "," << m_add_prior << ",";
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << -seine_temp3(i) << ",";
    }
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << -spawn_temp3(i) << ",";
    }
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << mdm_llk_ind(i) << ",";
    }
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << egg_llk_ind(i) << ",";
    }
    for (int i=1; i<=nyr_tobefit; i++){
        indiv_LLikReport << ADFG_hydro_llk_ind(i) << ",";
    }
    for (int i=1; i<=(nyr_tobefit); i++){
        indiv_LLikReport << PWSSC_hydro_llk_ind(i) << ",";
    }
    for (int i=1; i<=(nyr_tobefit-1); i++){
        indiv_LLikReport << juvenile_llk_ind(i) << ",";
    }
    indiv_LLikReport << juvenile_llk_ind(nyr_tobefit) << endl;
    

    // Disease outputs: Disease survival probs, infection indicences, fatality rates
    for (int i=1; i<=nyr_tobefit-1; i++){
      for (int j=1; j<=nage; j++){
        vhs_report1 << vhsv_survival_age(i,j) << ",";
      }

      vhs_report2 << inf_inc_sp(i) << ",";
      vhs_report3 << vhsprev_sp(i) << ",";
      vhs_report4 << fatal_sp(i) << ",";
      vhs_report5 << inf_inc_age3(i) << ",";
      vhs_report6 << vhsprev_age3(i) << ",";
      vhs_report7 << fatal_age3(i) << ",";

      for (int j=1; j<=nage; j++){
        ich_report1 << ich_survival_age(i,j) << ",";
      }
      ich_report2 << inf_inc_sp_ich(i) << ",";
      ich_report3 << ichprev_sp(i) << ",";
      ich_report4 << fatal_sp_ich(i) << ",";
      ich_report5 << inf_inc_age3_ich(i) << ",";
      ich_report6 << ichprev_age3(i) << ",";
      ich_report7 << fatal_age3_ich(i) << ",";
    }

    for (int j=1; j<=nage-1; j++){
        vhs_report1 << vhsv_survival_age(nyr_tobefit,j) << ","; 
    }
    vhs_report1 << vhsv_survival_age(nyr_tobefit,nage) << endl;
    vhs_report2 << inf_inc_sp(nyr_tobefit) << endl;
    vhs_report3 << vhsprev_sp(nyr_tobefit) << endl;
    vhs_report4 << fatal_sp(nyr_tobefit) << endl;
    vhs_report5 << inf_inc_age3(nyr_tobefit) << endl;
    vhs_report6 << vhsprev_age3(nyr_tobefit) << endl;
    vhs_report7 << fatal_age3(nyr_tobefit) << endl;
    
    for (int j=1; j<=nage-1; j++){
        ich_report1 << ich_survival_age(nyr_tobefit,j) << ","; 
    }
    ich_report1 << ich_survival_age(nyr_tobefit,nage) << endl;
    ich_report2 << inf_inc_sp_ich(nyr_tobefit) << endl;
    ich_report3 << ichprev_sp(nyr_tobefit) << endl;
    ich_report4 << fatal_sp_ich(nyr_tobefit) << endl;
    ich_report5 << inf_inc_age3_ich(nyr_tobefit) << endl;
    ich_report6 << ichprev_age3(nyr_tobefit) << endl;
    ich_report7 << fatal_age3_ich(nyr_tobefit) << endl;


    // Output PFRunBiomass for each saved iteration to .csv
    for (int i=1; i<=nyr_tobefit-1; i++){
     PFRReport << prefishery_biomass (i) << ","; 
    }
     PFRReport << prefishery_biomass (nyr_tobefit)<< ",";
     PFRReport << projected_prefishery_biomass << endl; // Projected pre-fishery run biomass for the upcoming year

    if(DD_Mat==1){
    	ofstream maturity_report("density_dependent_maturity.csv",ios::app);
    	for (int i=1; i<=nyr_tobefit; i++){
	      for (int j=1; j<=nage; j++){
	        maturity_report << maturity(i,j) << ",";
	      }
    	}
	}


REPORT_SECTION
  // report << "foo= " << endl << setprecision(10) << foo << endl; // to set precision just for foo OR
  // report.precision(10); // in first line of report section and to set precision of all report output
  // In order to get labels at the top of each of the following .csv I open them down here, write it once, 
  // then use ::app inside the mceval phase ios

  // Remove default report file output by ADMB because I want to contain everything in the rep_out folder
  remove("PWS_ASA.rep");

  struct stat buf;
  if(stat("rep_out",&buf)!=0){
    system("mkdir rep_out");
  }

  // User defined files
  ofstream SeACreport("rep_out/SeAC_pd.rep",ios::trunc);
  SeACreport << seine_age_comp << endl;
  SeACreport.close();
  
  ofstream SpACreport("rep_out/SpAC_pd.rep",ios::trunc);
  SpACreport << spawning_age_comp << endl;
  SpACreport.close();

  ofstream vhscompreport("rep_out/vhscomp_pd.rep",ios::trunc);
  vhscompreport << vhsv_pred << endl;
  vhscompreport.close();

  ofstream ichcompreport("rep_out/ichcomp_pd.rep",ios::trunc);
  ichcompreport << ich_pred << endl;
  ichcompreport.close();

  ofstream report1("rep_out/PWS_ASA.rep",ios::trunc);

  // the output below is sloppy - I need to display output for pre-1992 mort vectors without the disease index since that part of the vectors are place holders
  report1<<"LOG-LIKELIHOOD COMPONENTS" << endl;
  report1<< "penalty Count" << endl << penCount<< endl;
  report1<<"LL OF" << endl << full_llk << endl;
  report1<<"seine_llk" << endl << seine_llk <<endl;
  report1<<"spawner_llk" << endl << spawner_llk <<endl;
  report1<<"eggdep_llk " << endl << eggdep_llk <<endl;
  report1<<"ADFG_hydro_llk " << endl << ADFG_hydro_llk <<endl;
  report1<<"PWSSC_hydro_llk " << endl << PWSSC_hydro_llk <<endl;
  report1<<"mdm_llk " << endl <<mdm_llk<<endl;
  report1<<"age0_devs_penllk " << endl << age0_devs_penllk <<endl;
  report1<<"mort_devs_penllk " << endl << mort_devs_penllk <<endl;
  report1<<"age0_covar_prior " << endl << age0_covar_prior <<endl;
  report1<<"mort_covar_prior " << endl << mort_covar_prior <<endl;

  // Aerial juvenile survey (incorporated 12/2019)
  report1<<"juvenile_llk " << endl << juvenile_llk <<endl;

  // VHSV seroprevalence (11/2020)
  report1<<"vhsv_llk " << endl << vhsv_llk <<endl;

  // Ich. hoferi infection prevalence (11/2021)
  report1<<"ich_llk " << endl << ich_llk <<endl << endl;

  // Priors
  report1<<"Z_prior " << endl << prior_pars_1(3) <<endl;
  report1<<"mat_age3_prior " << endl << prior_pars_1(5) << endl;
  report1<<"mat_age4_prior " << endl << prior_pars_1(6) << endl;
  report1<<"hydADFG_add_prior " << endl << prior_pars_1(17) <<endl;
  report1<<"hydPWSSC_add_prior " << endl << prior_pars_1(19) <<endl;
  report1<<"m_add_prior " << endl << prior_pars_1(15) << endl;
  report1<<"vhsv_infection_prob prior " << endl << prior_pars_2(10) << endl;
  report1<<"ich_infection_prob prior " << endl << prior_pars_2(12) << endl;
  report1<<"vhsv_recovery_prob prior " << endl << prior_pars_2(11) << endl;
  report1<<"ich_recovery_prob prior " << endl << prior_pars_2(13) << endl << endl;

  report1<<"RESIDUALS" << endl;
  report1<<"Seine comps residuals" << endl << seine_comp_residuals << endl;
  report1<<"Spawner comps residuals" << endl << spawner_comp_residuals << endl;
  report1<<"Mile-days milt residuals" << endl << mdm_residuals << endl;
  report1<<"Egg deposition residuals" << endl << eggdep_residuals << endl;
  report1<<"ADFG Hydroacoustic residuals" << endl << ADFG_hydro_residuals << endl;
  report1<<"PWSSC Hydroacoustic residuals" << endl << PWSSC_hydro_residuals << endl << endl;

  report1<<"ANALYTICAL SIGMAS" << endl;
  report1<<"Combined Egg SD (egg_sd)" << endl <<egg_sd <<endl;
  report1<<"(Annual seine residuals)X(ESS)" << endl << seine_temp3 <<endl;
  report1<<"(Annual spawner residuals)X(ESS)" << endl << spawn_temp3 << endl << endl;

  
  report1 << "DERIVED QUANTITIES" << endl;
  report1 << "Pre-Fishery Run Biomass in mt" << endl << prefishery_biomass  << endl;
  report1 << "Pre-Fishery Run Biomass in TONS" << endl << prefishery_biomass_TONS  << endl;
  report1 << "Post-Fishery Spawning Biomass" << endl << postfishery_spawn_biomass << endl;
  report1 << "Estimated ADFG Hydro-acoustic Biomass" << endl << ADFG_HYDRO_est << endl;
  report1 << "Estimated PWSSC Hydro-acoustic Biomass" << endl << PWSSC_HYDRO_est << endl;
  report1 << "Pre-fishery total abundance (N_y_a)" << endl << N_y_a << endl;
  report1 << "Number of spawners (N_spawners_age)" << endl << N_spawners_age << endl << endl;

  report1 << "RECRUITMENT" << endl;
  report1 << "Recruits age-3" << endl;
  for (int i=1; i<=nyr_tobefit; i++){
    report1 << N_y_a(i,4) << endl; 
  }
  report1 << endl;

  report1 << "MATURITY" << endl;
  report1 << "Maturity-at-age of observed schools" << endl << maturity << endl << endl;
  report1 << "Maturity-at-age of unobserved schools" << endl << maturity_unobs << endl << endl;

  report1 << "ADULT SURVIVAL OUTPUTS" << endl;
  report1 << "Adult summer survival (summer_survival)" << endl << summer_survival << endl;
  report1 << "Adult winter survival (winter_survival)" << endl << winter_survival << endl << endl;
 
  report1 << "COVARIATE EFFECTS ON SURVIVAL" << endl;
  report1 << "Summer mortatlity" << endl << summer_mortality_effect << endl;
  report1 << "Winter mortality" << endl << winter_mortality_effect << endl << endl;

  report1 << "SUMMED ANNUAL MORTALITY DEVIATES (NON-ZERO IF ESTIMATED)" << endl;
  report1 << colsum(annual_mortdevs) << endl << endl;

  report1 << "ANNUAL MORTALITY DEVIATES (NON-ZERO IF ESTIMATED) - A MATRIX" << endl;
  report1 << annual_mortdevs << endl << endl;

  report1 << "VHSV ASSOCIATED SURVIVAL" << endl;
  report1 << "Proportion that avoids or survives infection:" << endl << vhsv_survival_age << endl << endl;

  report1 << "VHSV AGE-SPECIFIC IMMUNITY" << endl;
  report1 << "Proportions that are immune to reinfection:" << endl << vhsv_immunity_age << endl << endl;

  report1 << "VHSV OUTBREAK CHARACTERISTICS" << endl;
  report1 << "Incidence rate of infection in spawning population (numbers):" << endl << inf_inc_sp << endl << endl;
  report1 << "Fatality rate of infection in spawning population (numbers):" << endl << fatal_sp << endl << endl;
  report1 << "VHSV seroprevalence in spawning population (numbers):" << endl << vhsprev_sp << endl << endl;


  report1 << "ICH. HOF. INFECTION ASSOCIATED SURVIVAL" << endl;
  report1 << "Proportion that avoids or survives infection:" << endl << ich_survival_age << endl << endl;

  report1 << "ICH. HOF. INFECTION AGE-SPECIFIC CHRONIC INFECTION" << endl;
  report1 << "Proportions that are infected at start of year:" << endl << ich_chronic_infection_age << endl << endl;

  report1 << "ICH. HOF. INFECTION CHARACTERISTICS" << endl;
  report1 << "Incidence rate of infection in spawning population (numbers):" << endl << inf_inc_sp_ich << endl << endl;
  report1 << "Fatality rate of infection in spawning population (numbers):" << endl << fatal_sp_ich << endl << endl;
  report1 << "I. hoferi infection prevalence in spawning population (numbers):" << endl << ichprev_sp << endl << endl;

  //report1 << "PROJECTED MANAGEMENT QUANTITIES" << endl;
  //report1 << "Mean recruits from past 10 years" << endl << meanLgRec << endl;
  //report1 << "Projected total pre-fishery biomass" << endl << projected_prefishery_biomass << endl;
  //report1 << "Projected pre-fishery biomass by age" << endl << projected_Early_Sp_biomass << endl;
  //report1 << "Projected numbers at age" << endl << projected_N_y_a << endl << endl;

  report1.close();

RUNTIME_SECTION
maximum_function_evaluations 100, 100, 1000, 1000, 10000
