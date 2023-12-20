//compgV1 - Comprehensive Growth Model Version 1
//Data: combination of length frequency and mark-recapture data
//Model: estimates von Bertalanffy growth model parameters using a combination of length-frequency analysis and a growth increment model. Seasonal value for growth rate (k) are estiamted assuming growth is a linear function of growing degree days.
//Dec 20, 2023 
//Authors: Genny Nesslage, Mike Wilberg
//Acknowledgements: Some code for length-frequency analysis was borrowed from Dr. Robert Ahrens. Some tagging model code borrowed from Drs. Jie Cao, Yong Chen, and Larry Jacobson (file: lobster6f6.tpl 2013 American lobster assmt)
//Indexing notation: ts=time steps (nyrs*nseas),j=age,k=length bin

DATA_SECTION
	//.DAT file read in 
	init_int nseas //number of seasons
	init_int nsex //number of sexes
	init_int nage //guess for #ages
	init_int ts //number of time steps (nyrs*nseas)
	init_int fage //first age 
	init_int nbins //number of length bins
	init_int minbin //lower bound of smallest bin
	init_number binw //bin width (upper bound -lower bound)
	init_vector midbins(1,nbins) //middle of length bins
	//Environmental data for time-varyin K estimation
	init_vector GDD(1,ts) //vector of growing degree days
	
	//Length Frequency Analysis (LFA) submodel inputs
	init_int fEffN // effective sample size for females in LFA
	init_int mEffN // effective sample size for males in LFA
	init_matrix flfdata(1,ts,1,nbins) //length frequency data for females
	init_matrix mlfdata(1,ts,1,nbins) //length frequency data for males
	
	//Growth increment submodel inputs
	init_int   fTaggingDim  //Female tagging data dimension 
    init_int   mTaggingDim  //Male tagging data dimension
    init_matrix   fTaggingData(1,fTaggingDim,1,5) //Female tag data
    init_matrix   mTaggingData(1,mTaggingDim,1,5) //Male tag data  
	
	init_int eof //read-in error check
	!!if(eof!=999) cout<<"Error reading data file.\n Fix it."<<endl;
	!!if(eof!=999) ad_exit(1);
	
	//CONTROL FILE - switch to reading in starting values etc from .ctl file
	!!ad_comm::change_datafile_name("compgV1.ctl");
	
	//Common parameters - both submodels
    init_vector    fkpcprior(1,7)		 //growth parameter k proportionality constant (kpc; k is a linear function of kpc*GDD) for females
    init_vector    mkpcprior(1,7)		 //growth parameter k proportionality constant (kpc; k is a linear function of kpc*GDD) for males
	init_vector    fLinfprior(1,7)		 //growth parameter Linf for females
    init_vector    mLinfprior(1,7)		 //growth parameter Linf for males

	//LFA submodel parameters
	init_int lfirstph //phase for mean length for first age parameter
	init_number iflfirst // initial guess for mean length for first age female
	init_number imlfirst // initial guess for mean length for first age male
	init_matrix fsdl(1,ts,1,nage) // specified seasonal sd of mean length at age for females
	init_matrix msdl(1,ts,1,nage) // specified seasonal sd of mean length at age for males
	
	//Growth increment submodel parameters
    init_vector    fLinfSeprior(1,7)	 //sd of growth parameter Linf for females
    init_vector    mLinfSeprior(1,7)	 //sd of growth parameter Linf sd for males
    init_vector    fkSeprior(1,7)		 //sd of growth parameter k for females
    init_vector    mkSeprior(1,7)		 //sd of growth parameter k for males
    init_vector    fRhoprior(1,7)		 //correlation among growth parameters k and Linf for females
    init_vector    mRhoprior(1,7)		 //correlation among growth parameters k and Linf for males
	
	init_int eof2
	!!if(eof2!=999) cout<<"Error reading control file.\n Fix it."<<endl;
	!!if(eof2!=999)	ad_exit(1);

	vector age(1,nage); //vector of ages;
	!!age.fill_seqadd(fage,1.);

	int 		 phase     		//estimation phase
    number       lbd            //lower boundary
    number       ubd            //upper boundary
	int	j						//age index
	int k 						//bin index
	int t						//time step index
	
PARAMETER_SECTION
	//Set phase, lbd, upd, starting values for growth parameters common to BOTH SUBMODS by sex	  
	!!phase= (int) fkpcprior(1);lbd=fkpcprior(2);ubd=fkpcprior(3);
    init_bounded_number    log_fkpc(lbd,ubd,phase);
    !!phase= (int) mkpcprior(1);lbd=mkpcprior(2);ubd=mkpcprior(3);
    init_bounded_number    log_mkpc(lbd,ubd,phase);
	!!phase= (int)fLinfprior(1);lbd=fLinfprior(2);ubd=fLinfprior(3);
    init_bounded_number    log_fLinf(lbd,ubd,phase);
    !!phase= (int) mLinfprior(1);lbd=mLinfprior(2);ubd=mLinfprior(3);
    init_bounded_number    log_mLinf(lbd,ubd,phase);
	!!phase= (int) fkSeprior(1);lbd=fkSeprior(2);ubd=fkSeprior(3);
    init_bounded_number    log_fkSe(lbd,ubd,phase);
    !!phase= (int) mkSeprior(1);lbd=mkSeprior(2);ubd=mkSeprior(3);
    init_bounded_number    log_mkSe(lbd,ubd,phase);
    !!phase= (int) fLinfSeprior(1);lbd=fLinfSeprior(2);ubd=fLinfSeprior(3);
    init_bounded_number    log_fLinfSe(lbd,ubd,phase);
    !!phase= (int) mLinfSeprior(1);lbd=mLinfSeprior(2);ubd=mLinfSeprior(3);
    init_bounded_number    log_mLinfSe(lbd,ubd,phase);
    !!phase= (int) fRhoprior(1);lbd=fRhoprior(2);ubd=fRhoprior(3);
    init_bounded_number    fRho(lbd,ubd,phase);
    !!phase= (int) mRhoprior(1);lbd=mRhoprior(2);ubd=mRhoprior(3);
    init_bounded_number    mRho(lbd,ubd,phase);
	
	//Estimate growth parameters in log space
	!!log_fkpc=log(fkpcprior(7));
	!!log_mkpc=log(mkpcprior(7));
	!!log_fLinf=log(fLinfprior(7));
	!!log_mLinf=log(mLinfprior(7));
	!!log_fkSe=log(fkSeprior(7));
	!!log_mkSe=log(mkSeprior(7));
	!!log_fLinfSe=log(fLinfSeprior(7));
	!!log_mLinfSe=log(mLinfSeprior(7));
		
	//LFA-specific parameters
	init_number log_flfirst(lfirstph); //mean length for first age in first season female
	init_number log_mlfirst(lfirstph); //mean length for first age in first season male
	init_bounded_matrix log_fnpage(1,ts,fage+1,nage,-10.,10.); //numbers at age for females ages 2+
	init_bounded_matrix log_mnpage(1,ts,fage+1,nage,-10.,10.); //numbers at age for males ages 2+

	//Estimate LFA parameters in log space
	!!log_flfirst=log(iflfirst);
	!!log_mlfirst=log(imlfirst);
			
	//Storage for estimated and calculated values
	vector fk(1,ts) //female k in each time step
	vector mk(1,ts) //male k in each time step
	number fkSe //female sd of k
	number mkSe //male sd of k
	number fLinf //female Linf (asymptotic length)
	number mLinf //male Linf (asymptotic length)
	number fLinfSe //female sd of Linf
	number mLinfSe //male sd of Linf	
	
	//LFA submodel
	number ffyrk //female k in the first year (sum across first nseas)
	number mfyrk //male k in the first year (sum across first nseas)
	number flfirst //mean length of first age female
	number mlfirst //mean length of first age male
	matrix fnpage(1,ts,1,nage) //female numbers at age
	matrix mnpage(1,ts,1,nage) //male numbers at age
	vector fpl(1,nbins); //female proportions at length (LFA)
	vector mpl(1,nbins); //male proportions at length
	matrix fmeanl(1,ts,1,nage) //female mean length at age
	matrix mmeanl(1,ts,1,nage) //male mean length at age
	matrix fpla(1,nbins,1,nage) //female probability of length at age
	matrix mpla(1,nbins,1,nage) //male probability of length at age
	matrix fppl(1,ts,1,nbins) //female proportions at length
	matrix mppl(1,ts,1,nbins) //male proportions at length

	//storage for Growth Increment (GI) submodel
	vector  fDeltaL(1,fTaggingDim) //observed female growth between capture and recapture
	vector  mDeltaL(1,mTaggingDim) //observed male growth between capture and recapture
	vector  fDeltaT(1,fTaggingDim) //observed female time at large
	vector  mDeltaT(1,mTaggingDim) //observed male time at large
	vector  fLt1(1,fTaggingDim)    //observed female length at capture
	vector  mLt1(1,mTaggingDim)    //observed male length at capture
	vector  ftcap(1,fTaggingDim)   //observed female time step at capture
	vector  mtcap(1,mTaggingDim)   //observed male time step at capture
	vector  ftrecap(1,fTaggingDim) //observed female time step at recapture 
	vector  mtrecap(1,mTaggingDim) //observed male time step at recapture
	vector  fDeltaLhat(1,fTaggingDim) //store female predicted growth increments
	vector  mDeltaLhat(1,mTaggingDim) //store male predicted growth increments
    number fDeltaSD; //standard deviation for female tagging data residuals
	number mDeltaSD; //standard deviation for male tagging data residuals
	vector  fDeltaStdResid(1,fTaggingDim) //female standardized residuals for tagging estimates
	vector  mDeltaStdResid(1,mTaggingDim) //male standardized residuals for tagging estimates
	vector  ftagsd(1,fTaggingDim);  //store female sd of growth increments
	vector  mtagsd(1,mTaggingDim);  //store male sd of growth increments
 	 
	//storage for all seasonal growth transition matrices (GTMs)
	3darray FG(1,ts,1,nbins,1,nbins)
	3darray MG(1,ts,1,nbins,1,nbins)

	number comp1
	number comp2
	number comp3
	number comp4
	//number npar
	//number AIC
	//number AICc
	//number BIC	

	objective_function_value nll;
	number grad_max;

PRELIMINARY_CALCS_SECTION
	if (fTaggingDim>0)
     {
		 fDeltaL=column(fTaggingData,5);
		 fDeltaT=column(fTaggingData,4);
		 fLt1=column(fTaggingData,3);
		 ftcap=column(fTaggingData,1);
		 ftrecap=column(fTaggingData,2);
	 }

	 if (mTaggingDim>0)
	 {
		 mDeltaL=column(mTaggingData,5);
		 mDeltaT=column(mTaggingData,4);
		 mLt1=column(mTaggingData,3);
		 mtcap=column(mTaggingData,1);
		 mtrecap=column(mTaggingData,2); 
	 }

	for(int t=1;t<=ts;t++){
		log_fnpage(t)=0.0;
		log_mnpage(t)=0.0;
	}	  

	 
PROCEDURE_SECTION
	//initialize parameters
	initialization();
	
	//generate sex-specific GTMs for each time step
	for(int t=1;t<=ts;t++)
	{
		FG(t)=CalGrowthMatrixSimple(fLinf,fk(t),fLinfSe,fkSe,fRho);
		MG(t)=CalGrowthMatrixSimple(mLinf,mk(t),mLinfSe,mkSe,mRho);	
	}
	
	//LFA submodel
	lfa_calculations();

	//Growth increment submodel for each sex
	fGI_model();
	mGI_model();
	
	//Residuals and standardized residuals for tagging data
	//calc GI model residuals for objective function
		fDeltaStdResid=fDeltaL-fDeltaLhat;
        mDeltaStdResid=mDeltaL-mDeltaLhat; 
	//standard deviations for residuals
        fDeltaSD=sqrt(norm2(fDeltaStdResid)/double(fTaggingDim));
        mDeltaSD=sqrt(norm2(mDeltaStdResid)/double(mTaggingDim));
	//standardized residuals
        fDeltaStdResid/=fDeltaSD;
        mDeltaStdResid/=mDeltaSD;
			
	objective_function();
	
	if(last_phase())
	{
		grad_max=objective_function_value::pobjfun->gmax;
	} 		
	
FUNCTION initialization
	//Common parameters
	fk=GDD*mfexp(log_fkpc);
	mk=GDD*mfexp(log_mkpc);
	fLinf=mfexp(log_fLinf);
	mLinf=mfexp(log_mLinf);
	fkSe=mfexp(log_fkSe);
	mkSe=mfexp(log_mkSe);
	fLinfSe=mfexp(log_fLinfSe);
	mLinfSe=mfexp(log_mLinfSe);
	fRho=fRhoprior(7);
	mRho=mRhoprior(7);

	//LFA parameters
	ffyrk=fk[1];
	mfyrk=fk[1];
	for(t=2;t<=nseas;t++) 
		{
			ffyrk+=fk[t];
			mfyrk+=mk[t];
		}	
	flfirst=mfexp(log_flfirst);
	mlfirst=mfexp(log_mlfirst);

	
FUNCTION dvar_matrix CalGrowthMatrixSimple(dvariable Linf,prevariable K,dvariable LSigma,dvariable KSigma,dvariable Rho)
     
	// Arguments:  Linf=maximum size (L-infinity)
	//             K=von Bertlanffy growth parameter **vector**
	//             LSigma=standard deviation for Linf in population
	//             KSigma=standard deviation for K in population
	//             Rho=correlation between Linf and K in popn (constrained at pres as -1 < rho < 1)
	//			   Cov = Rho*LSigma*KSigma
	// Output: square growth transition matrix with starting sizes in rows and ending sizes
	//         in columns by the method of Chen et al. (2003)
	// Originally coded by Jie Cao (Univ. of Maine) and revised by Larry Jacobson (Nov. 19, 2013) and Genny Nesslage/Mike Wilberg (2023)

		// place to store function output
		dvar_matrix TransitionMatrix(1,nbins,1,nbins);
		TransitionMatrix.initialize();

		// used repeatedly later
		dvariable eKbar=mfexp(-K);
		dvariable meKbar=1.0-eKbar;

		// mean and variance for growth increments
		dvariable deltaLbar;
		dvariable SddeltaLbar;

		// use middle of bins to characterize starting length
		// n is the starting bin
			for (int n=1; n<nbins; n++) 
				{
				// middle of starting bin relative to Linf
				dvariable dmidsize=Linf-midbins(n);
			
				// mean increment at this starting size
					deltaLbar=dmidsize*meKbar;
				// mean size after growth for this starting size
					dvariable newMeanLbar=midbins(n)+deltaLbar;

				// variance of increments Eqn 14 in Cao et al. 2016a
					SddeltaLbar=square(LSigma*meKbar); //first term
					SddeltaLbar+=KSigma*square(dmidsize*eKbar); // middle term
					SddeltaLbar+=2.0*Rho*KSigma*LSigma*meKbar*dmidsize*eKbar; //last term; -= constrains covariance of Linf and k to be neg
					SddeltaLbar=sqrt(SddeltaLbar); //calc sd from variance
		
					for (int d=n; d<nbins; d++)
					{
			
					// lower and upper bounds for receiving bin
					// standardize to z-scores	
						  dvariable lb=((midbins(d)-0.5*binw)-newMeanLbar)/SddeltaLbar;
						  dvariable ub=((midbins(d+1)-0.5*binw)-newMeanLbar)/SddeltaLbar;

					// probability
						  TransitionMatrix(n,d)=cumd_norm(ub)-cumd_norm(lb);
						}
		
				// last size bin and normalization
					TransitionMatrix(n)/=sum(TransitionMatrix(n));
				//  if (n>=nbins) exit(1);
				  }
			//EO(TransitionMatrix);
			//exit(0);
	  return(TransitionMatrix);
		  	  
FUNCTION lfa_calculations
				
		//Mean LAA for first age in first season
		
		fmeanl(1,1)=flfirst;
		mmeanl(1,1)=mlfirst;
					
		//Mean LAA for all other ages in first season 
		
		for(int j=2;j<=nage;j++) 
			{
				fmeanl(1,j)=fmeanl(1,j-1)+(fLinf-fmeanl(1,j-1))*(1.-exp(-ffyrk));			
				mmeanl(1,j)=mmeanl(1,j-1)+(mLinf-mmeanl(1,j-1))*(1.-exp(-mfyrk));
			}
		
		//Mean LAA for all ages in laters seasons
		for(int t=2;t<=ts;t++)
		{            
			if(t%nseas==1)
			{
				//first age mean length
				fmeanl(t,1)=flfirst;
				mmeanl(t,1)=mlfirst;
				//rest of the ages
				for(int j=2;j<=nage;j++)
				{ 
					fmeanl(t,j)=fmeanl(t-1,j-1)+(fLinf-fmeanl(t-1,j-1))*(1.-exp(-fk(t-1)));
					mmeanl(t,j)=mmeanl(t-1,j-1)+(mLinf-mmeanl(t-1,j-1))*(1.-exp(-mk(t-1)));
				}
			}
			else
			{		
				for(int j=1;j<=nage;j++)
				{ 	
					fmeanl(t,j)=fmeanl(t-1,j)+(fLinf-fmeanl(t-1,j))*(1.-exp(-fk(t-1)));
					mmeanl(t,j)=mmeanl(t-1,j)+(mLinf-mmeanl(t-1,j))*(1.-exp(-mk(t-1)));
				}
			}
		}

		//Calculate numbers at age
		for(int t=1;t<=ts;t++)
			{  
			fnpage(t,1)=1.; //props at age est'd relative to first 
			mnpage(t,1)=1.;
				for(int j=2;j<=nage;j++)
				{
					fnpage(t,j)=mfexp(log_fnpage(t,j));
					mnpage(t,j)=mfexp(log_mnpage(t,j));
				}
				fnpage(t)=fnpage(t)/sum(fnpage(t));
				mnpage(t)=mnpage(t)/sum(mnpage(t));
			}
	
		dvariable fz1;
		dvariable fz2;
		dvariable mz1;	
		dvariable mz2;	
		
		//calc length probs and props
		for(int t=1;t<=ts;t++) 
			{
				for(int j=1;j<=nage;j++) //loop over ages to calc length probs
				{
					for(int k=1;k<=nbins;k++) //loop over length bins
					{
						fz1=((midbins(k)-0.5*binw)-fmeanl(t,j))/fsdl(t,j); 
						mz1=((midbins(k)-0.5*binw)-mmeanl(t,j))/msdl(t,j); 
						fz2=((midbins(k)+0.5*binw)-fmeanl(t,j))/fsdl(t,j);
						mz2=((midbins(k)+0.5*binw)-mmeanl(t,j))/msdl(t,j);
						fpla(k,j)=cumd_norm(fz2)-cumd_norm(fz1); 		
						mpla(k,j)=cumd_norm(mz2)-cumd_norm(mz1); 
					}
				}
			//estimate proportions in each length bin
			fpl=fpla*fnpage(t); 
			mpl=mpla*mnpage(t); 
			fppl(t)=fpl/sum(fpl);  //standardize
			mppl(t)=mpl/sum(mpl);  //standardize
			}


FUNCTION dvariable varGI(dvariable length,dvariable Linf,dvariable K, dvariable LSigma,dvariable KSigma, dvariable Rho) 
		dvariable eKbar=mfexp(-K);
		dvariable meKbar=1.0-eKbar;
		dvariable vardeltaLbar;
		dvariable dmidsize=Linf-length;
	
		vardeltaLbar=square(LSigma*meKbar); //first term
		vardeltaLbar+=KSigma*square(dmidsize*eKbar); // middle term
		vardeltaLbar+=2.0*Rho*LSigma*KSigma*meKbar*dmidsize*eKbar; 
		return(vardeltaLbar); 
		

FUNCTION fGI_model
	for(int r=1;r<=fTaggingDim;r++){
		dvariable tmpDeltaLhat=0.0;	
		dvariable tmpDeltaLhatstore=0.0;	
		dvariable tmpL=0.0;
		dvariable tmpvar=0.0;
		
		int capts=int(fTaggingData(r,1));
		int recapts=int(fTaggingData(r,2));
		tmpL=fTaggingData(r,3);

		for(int t=capts;t<recapts;t++){
			if(t==capts){
				tmpDeltaLhat=(fLinf-tmpL)*(1.0-mfexp(-fk(t)*1.0));
				tmpDeltaLhatstore=tmpDeltaLhat;
				tmpvar=varGI(tmpL,fLinf,fk(t),fLinfSe,fkSe,fRho);
						}
			else{
				tmpL+=tmpDeltaLhat;
				tmpDeltaLhat=(fLinf-tmpL)*(1.0-mfexp(-fk(t)*1.0));
				tmpDeltaLhatstore+=tmpDeltaLhat;
				tmpvar+=varGI(tmpL,fLinf,fk(t),fLinfSe,fkSe,fRho);			
				}
			}
		fDeltaLhat(r)=tmpDeltaLhatstore;
		ftagsd(r)=sqrt(tmpvar);
		}

						
FUNCTION mGI_model
	for(int r=1;r<=mTaggingDim;r++){
		dvariable tmpDeltaLhat=0.0;	
		dvariable tmpDeltaLhatstore=0.0;
		dvariable tmpL=0.0;
		dvariable tmpvar=0.0;

		int capts=int(mTaggingData(r,1));
		int recapts=int(mTaggingData(r,2));
		tmpL=mTaggingData(r,3);
		
		for(int t=capts;t<recapts;t++){
			if(t==capts){
				tmpDeltaLhat=(mLinf-tmpL)*(1.0-mfexp(-mk(t)*1.0));
				tmpDeltaLhatstore=tmpDeltaLhat;
				tmpvar=varGI(tmpL,mLinf,mk(t),mLinfSe,mkSe,mRho);
					}
			else{
				tmpL+=tmpDeltaLhat;
				tmpDeltaLhat=(mLinf-tmpL)*(1.0-mfexp(-mk(t)*1.0));
				tmpDeltaLhatstore+=tmpDeltaLhat;
				tmpvar+=varGI(tmpL,mLinf,mk(t),mLinfSe,mkSe,mRho);	
				}
				//cout<<tmpvar<<endl;
			}      		
		mDeltaLhat(r)=tmpDeltaLhatstore;
		mtagsd(r)=sqrt(tmpvar);
		}	
					
FUNCTION objective_function 

	//Multinomial
	//LFA components
	comp1=0.;
	comp1=-double(fEffN)*sum(elem_prod(flfdata,log(fppl+0.0001))); 
	comp2=-double(mEffN)*sum(elem_prod(mlfdata,log(mppl+0.0001))); 
	// likelihood for growth increment model - normal distribution
	comp3=sum(log(ftagsd))+0.5*norm2(elem_div((fDeltaL-fDeltaLhat),ftagsd));   
	comp4=sum(log(mtagsd))+0.5*norm2(elem_div((mDeltaL-mDeltaLhat),mtagsd)); 

	nll=comp1+comp2+comp3+comp4;  
	
	//total objective function
	//npar=nseas*nage+6.;
	//AIC=2*npar-2.*(comp1+comp2+comp3);
	//AICc=AIC+(2*npar*(npar-1.))/(nseas*double(nbins)-npar-1.);
	//BIC=npar*log(nseas*double(nbins))-2.*(comp1+comp2+comp3);

REPORT_SECTION
  report<<"femaleK "<<fk<<endl;
  report<<"maleK "<<mk<<endl;
  report<<"femaleLinf "<<fLinf<<endl;
  report<<"maleLinf "<<mLinf<<endl;
  report<<"female_mnL@fage "<<flfirst<<endl;
  report<<"male_mnL@fage "<<mlfirst<<endl;
  report<<"female_npage "<<endl;
  report<<fnpage<<endl;
  report<<"male_npage"<<endl;
  report<<mnpage<<endl;
 
  ofstream ofs1("LFA_est.rep");
	ofs1<<fppl<<endl;
	ofs1<<mppl<<endl;
  
  ofstream ofs2("LFA_obs.rep");
	ofs2<<flfdata<<endl;
	ofs2<<mlfdata<<endl;

  ofstream ofs3("GI_est.rep");
	ofs3<<fDeltaLhat<<endl;
	ofs3<<mDeltaLhat<<endl;
  
  ofstream ofs4("GI_obs.rep");
	ofs4<<fDeltaL<<endl;
	ofs4<<mDeltaL<<endl;
	
  ofstream ofs5("tvK.rep");
	ofs5<<fk<<endl;
	ofs5<<mk<<endl;
	
	
TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	#define EO(object) cout << #object "= " << object << endl << endl

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"grad_max"<<grad_max<<endl;
	cout<<"female vbk"<<fk<<endl;
	cout<<"male vbk"<<mk<<endl;
	cout<<"female linf"<<fLinf<<endl;
	cout<<"male linf"<<mLinf<<endl;
	cout<<"female lfirst"<<flfirst<<endl;
	cout<<"male lfirst"<<mlfirst<<endl;
	cout<<comp1<<endl;
	cout<<comp2<<endl;
	cout<<comp3<<endl;
	cout<<comp4<<endl;
	cout<<nll<<endl;
	cout<<"*******************************************"<<endl;


